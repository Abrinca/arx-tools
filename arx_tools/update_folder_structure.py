import os
import sys
import json
import logging
import shutil
import tarfile
import warnings

from .check_v3 import check_genome_v3
from .folder_looper import FolderLooper, FolderGenome
from .rename_eggnog import EggnogFile
from .rename_genbank import GenBankFile
from .rename_fasta import FastaFile
from .utils import query_yes_no, get_folder_structure_version


def _get_folder_structure_dir(folder_structure_dir: str = None) -> str:
    if folder_structure_dir is None:
        if 'FOLDER_STRUCTURE' not in os.environ:
            raise ValueError('Cannot find the folder_structure. Please set --folder_structure_dir or environment variable FOLDER_STRUCTURE')
        folder_structure_dir = os.environ['FOLDER_STRUCTURE']
    if not os.path.isdir(folder_structure_dir):
        raise ValueError(f'Could not find the folder_structure. Folder does not exist: {folder_structure_dir}')
    return folder_structure_dir


def set_folder_structure_version(new_version: int, folder_structure_dir: str) -> None:
    if not isinstance(folder_structure_dir, str):
        raise TypeError(f'folder_structure_dir must be str, got {type(folder_structure_dir).__name__}')
    version_file = f'{folder_structure_dir}/version.json'

    with open(version_file) as f:
        version_dict = json.load(f)

    version_dict['folder_structure_version'] = new_version

    with open(version_file, 'w') as f:
        json.dump(version_dict, f, indent=4)

    print()
    print(f'Successfully updated to folder structure version {new_version}!')


def ask(v_from: int, v_to: int, actions: [str], folder_structure_dir: str):
    if not isinstance(folder_structure_dir, str):
        raise TypeError(f'folder_structure_dir must be str, got {type(folder_structure_dir).__name__}')

    current_version = get_folder_structure_version(folder_structure_dir)
    if current_version != v_from:
        raise ValueError(
            f'Cannot proceed: Folder structure version mismatch.\n'
            f'This script expects version {v_from}, but folder_structure/version.json says version {current_version}.'
        )

    question = f'Upgrade folder structure from version {v_from} to {v_to}:'
    for action in actions:
        question += f'\n - {action}'
    question += '\n\nProceed?'
    if not query_yes_no(question=question, default=None):
        exit(1)


def loop_genomes(folder_structure_dir: str, skip_ignored=False, sanity_check=False, representatives_only=False) -> [FolderGenome]:
    for genome in FolderLooper(folder_structure_dir=folder_structure_dir).genomes(
            skip_ignored=skip_ignored,
            sanity_check=sanity_check,
            representatives_only=representatives_only
    ):
        try:
            has_json = genome.has_json
        except PermissionError:
            print(f'{genome.identifier}: SKIPPED: permission denied on genome.json')
            continue
        if has_json:
            yield genome


def from_1_to_2(folder_structure_dir: str = None, skip_ignored=False, sanity_check=False, representatives_only=False):
    """ Upgrade OpenGenomeBrowser folder structure. """
    folder_structure_dir = _get_folder_structure_dir(folder_structure_dir)
    v_from = 1
    v_to = 2

    ask(v_from=v_from, v_to=v_to, actions=['add COG to genome.json'], folder_structure_dir=folder_structure_dir)

    for genome in loop_genomes(folder_structure_dir=folder_structure_dir, skip_ignored=skip_ignored, sanity_check=sanity_check,
                               representatives_only=representatives_only):
        genome_json = genome.json
        if 'COG' in genome_json:
            print(f'{genome.identifier}: already has COG in genome.json')
            continue

        cog = {}  # default

        eggnog_files = [f for f in genome_json['custom_annotations'] if f['type'].startswith('eggnog')]
        for file in eggnog_files:
            path = os.path.join(genome.path, file['file'])
            try:
                cog = EggnogFile(file=path).cog_categories()
            except AssertionError as e:
                logging.info(msg=str(e))
                pass

        print(f'{genome.identifier}: adding COG={cog}')
        genome_json['COG'] = cog
        genome.replace_json(genome_json)

    set_folder_structure_version(new_version=v_to, folder_structure_dir=folder_structure_dir)



def _apply_contig_map_to_fna(src: str, dst: str, contig_map: dict) -> int:
    """Write a contig-header-renamed copy of a FASTA file to dst. Returns rename count."""
    with open(src) as f:
        lines = f.readlines()

    if not lines:
        raise ValueError(f"FASTA file is empty: {src}")

    renamed = 0
    result = []
    for line in lines:
        if line.startswith('>'):
            parts = line[1:].split(None, 1)
            old_id = parts[0]
            new_id = _resolve_contig_id(old_id, contig_map)
            if new_id is None:
                raise ValueError(
                    f'{os.path.basename(src)}: contig {old_id!r} not found in contig_map '
                    f'(GBK has {list(contig_map)[:3]}...); '
                    f'assembly FNA and GBK contig IDs may be out of sync')
            rest = (' ' + parts[1]) if len(parts) > 1 else '\n'
            line = f'>{new_id}{rest}'
            renamed += 1
        result.append(line)

    with open(dst, 'w') as f:
        f.writelines(result)
    return renamed


def _rename_sequence_region(line: str, contig_map: dict) -> str:
    """Rename the seqid in a ##sequence-region pragma line."""
    parts = line.split(None, 3)
    if len(parts) != 4:
        return line
    new_contig = _resolve_contig_id(parts[1], contig_map)
    if new_contig is None:
        return line
    parts[1] = new_contig
    line = ' '.join(parts)
    if not line.endswith('\n'):
        line += '\n'
    return line


def _rename_gff_feature_line(line: str, contig_map: dict) -> tuple[str, int]:
    """Rename seqid in a GFF3 feature line; return (line, seqid_changed)."""
    cols = line.split('\t')
    if len(cols) != 9:
        return line, 0
    new_contig = _resolve_contig_id(cols[0], contig_map)
    if new_contig is None:
        return line, 0
    cols[0] = new_contig
    return '\t'.join(cols), 1


def _apply_contig_map_to_gff(src: str, dst: str, contig_map: dict) -> int:
    """Write a contig-seqid-renamed copy of a GFF file to dst. Returns seqid_changed count."""
    with open(src) as f:
        lines = f.readlines()

    seqid_changed = 0
    in_fasta = False
    result = []
    for line in lines:
        if in_fasta:
            if line.startswith('>'):
                parts = line[1:].split(None, 1)
                new_id = _resolve_contig_id(parts[0], contig_map)
                if new_id is not None:
                    rest = (' ' + parts[1]) if len(parts) > 1 else '\n'
                    line = f'>{new_id}{rest}'
        elif line.strip() == '##FASTA':
            in_fasta = True
        elif line.startswith('##sequence-region'):
            line = _rename_sequence_region(line, contig_map)
        elif not line.startswith('#') and line.strip():
            line, sc = _rename_gff_feature_line(line, contig_map)
            seqid_changed += sc
        result.append(line)

    with open(dst, 'w') as f:
        f.writelines(result)
    return seqid_changed


def _strip_gnl_prefix(contig_id: str) -> str:
    """Strip gnl|X| or similar pipe-delimited prefix, returning the bare ID."""
    return contig_id.rsplit('|', 1)[1] if '|' in contig_id else contig_id


def _strip_version_suffix(contig_id: str) -> str:
    """Strip NCBI version suffix (.N) from a contig ID, returning the bare accession."""
    if '.' in contig_id:
        base, suffix = contig_id.rsplit('.', 1)
        if suffix.isdigit():
            return base
    return contig_id


def _resolve_contig_id(contig_id: str, contig_map: dict) -> str | None:
    """Resolve FNA contig ID against GBK-keyed map, tolerating gnl|X| prefixes and version mismatches."""
    new_id = contig_map.get(contig_id)
    if new_id is not None:
        return new_id
    base = _strip_gnl_prefix(contig_id)
    if base != contig_id and base in contig_map:
        return contig_map[base]
    for key, val in contig_map.items():
        if _strip_version_suffix(key) == contig_id:
            return val
    return None



def _locus_tags_need_rename(gbk_path: str, genome_id: str) -> bool:
    """
    Return True if the GBK contains locus tags but none start with '{genome_id}_'.

    Scans the file as text (no BioPython parse) so it's fast even for large GBKs.
    Returns False if locus tags are already arx-assigned, or if there are no locus tags.
    """
    prefix = f'/locus_tag="{genome_id}_'
    found_any = False
    with open(gbk_path) as f:
        for line in f:
            s = line.lstrip()
            if s.startswith('/locus_tag="'):
                found_any = True
                if s.startswith(prefix):
                    return False
    return found_any


_BLAST_EXTENSIONS = {
    # nucleotide v4/v5
    '.nhr', '.nin', '.nsq', '.nsi', '.nsd', '.ndb', '.not', '.ntf', '.nto', '.njs',
    # protein v4/v5
    '.phr', '.pin', '.psq', '.psi', '.psd', '.pdb', '.pot', '.ptf', '.pto', '.pjs',
    # alias files
    '.nal', '.pal',
}


def _promote_v3_files(v3_to_orig: dict, genome_dir: str, genome_id: str,
                      extra_backup: list = None) -> str:
    """Archive originals into {genome_id}_v2_backup.tar.gz and rename .v3 files into place; return archive path."""
    archive_path = os.path.join(genome_dir, f'{genome_id}_v2_backup.tar.gz')
    with tarfile.open(archive_path, 'w:gz') as tar:
        for orig_path in v3_to_orig.values():
            if os.path.exists(orig_path):
                tar.add(orig_path, arcname=os.path.relpath(orig_path, genome_dir))
        for extra_path in (extra_backup or []):
            if os.path.exists(extra_path):
                tar.add(extra_path, arcname=os.path.relpath(extra_path, genome_dir))

    for v3_path, orig_path in v3_to_orig.items():
        if os.path.exists(orig_path):
            os.unlink(orig_path)
        shutil.move(v3_path, orig_path)

    return archive_path


def from_2_to_3(folder_structure_dir: str = None, skip_ignored=False, contig_format: str = '_scf{n}',
                create_from_file: bool = False, create_only: bool = False, promote: bool = False):
    """
    Upgrade folder structure from v2 to v3.

    Only contig IDs are updated; locus tags are left untouched (arx enforces correct locus tags
    at import time, so all genomes in the folder structure already have the right prefix).

    Per genome:
      1. Shallow v3 check: skip if already v3; warn if a partial upgrade (.v3 files) exists.
      2. Generate .v3 intermediate files (gbk, assembly fna, gff).
         Originals are untouched during this phase.
      3. Archive originals into {genome_id}_v2_backup.tar.gz and promote .v3 files into place.
         With --create_from_file, regenerate gff/faa/ffn from the updated GBK instead.
      4. Post-check to verify success.
      5. Delete BLAST databases (they reference stale contig IDs; rebuild manually in arx).

    Two-step workflow:
      --create_only  Stop after step 2: generate .v3 files without promoting them.
                     Inspect the generated files, then re-run with --promote to finish.
      --promote      Skip step 2: promote .v3 files left by a previous --create_only run.
    """
    folder_structure_dir = _get_folder_structure_dir(folder_structure_dir)

    warnings.filterwarnings('ignore', message='.*malformed locus line.*', module='Bio')
    warnings.filterwarnings('ignore', message='.*Premature end of file.*', module='Bio')
    warnings.filterwarnings('ignore', message='.*Expected sequence length.*', module='Bio')

    if create_only:
        actions = [
            'shallow-check each genome; skip if already v3',
            'generate .v3 for gbk, assembly fna, and gff (contig IDs only)',
            '(promotion skipped; re-run with --promote to archive originals and promote)',
        ]
    elif promote:
        actions = [
            'shallow-check each genome; skip if already v3 or no pending .v3 files',
            'archive originals into {genome_id}_v2_backup.tar.gz and promote pending .v3 files',
            'post-check each genome to verify',
            'delete BLAST databases (rebuild manually in arx when needed)',
        ]
    elif create_from_file:
        actions = [
            'shallow-check each genome; skip if already v3',
            'generate .v3 for gbk and assembly fna (contig IDs only)',
            'archive originals and regenerate gff/faa/ffn from updated GBK',
            'post-check each genome to verify',
            'delete BLAST databases (rebuild manually in arx when needed)',
        ]
    else:
        actions = [
            'shallow-check each genome; skip if already v3',
            'generate .v3 for gbk, assembly fna, and gff (contig IDs only)',
            'archive originals into {genome_id}_v2_backup.tar.gz and promote .v3 files',
            'post-check each genome to verify',
            'delete BLAST databases (rebuild manually in arx when needed)',
        ]
    ask(v_from=2, v_to=3, actions=actions, folder_structure_dir=folder_structure_dir)

    succeeded = 0
    skipped_not_ready = 0   # --promote: genomes still v2 with no pending .v3 files
    skipped_pending = 0     # normal run: genomes with leftover .v3 files (warn and skip)
    skipped_needs_rename = 0  # genomes with external locus tags that need arx --rename first
    skipped_errors = 0        # genomes skipped due to missing/misconfigured files

    for genome in loop_genomes(folder_structure_dir=folder_structure_dir, skip_ignored=skip_ignored):
        genome_json = genome.json
        genome_id = genome.identifier

        # 1. Shallow check
        pre_check = check_genome_v3(genome.path, genome_id, deep=False, contig_format=contig_format)
        if pre_check.is_v3:
            print(f'{genome_id}: already v3, skipping')
            continue

        if promote:
            # --promote: only handle .v3 files left by a previous --create_only run.
            if not pre_check.has_pending_v3_files:
                print(f'{genome_id}: no pending .v3 files, skipping (run with --create_only first)')
                skipped_not_ready += 1
                continue
            names = ', '.join(os.path.basename(p) for p in pre_check.pending_files)
            print(f'{genome_id}: promoting {names}')
            v3_to_orig = {p: p[:-len('.v3')] for p in pre_check.pending_files}
            archive = _promote_v3_files(v3_to_orig, genome_dir=genome.path, genome_id=genome_id)
            print(f'{genome_id}: archived originals → {os.path.basename(archive)}')
            post_check = check_genome_v3(genome.path, genome_id, deep=False, contig_format=contig_format)
            if not post_check.is_v3:
                reasons = '; '.join(post_check.issues)
                print(f'\n{"=" * 60}', file=sys.stderr)
                print(f'ERROR: {genome_id}: post-check failed after promote', file=sys.stderr)
                for issue in post_check.issues:
                    print(f'  - {issue}', file=sys.stderr)
                print(f'{"=" * 60}\n', file=sys.stderr)
                raise ValueError(f'{genome_id}: post-check failed after promote: {reasons}')
            print(f'{genome_id}: done (post-check OK)')
            succeeded += 1
            deleted = 0
            for fname in os.listdir(genome.path):
                if any(fname.endswith(ext) for ext in _BLAST_EXTENSIONS):
                    os.remove(os.path.join(genome.path, fname))
                    deleted += 1
            if deleted:
                print(f'{genome_id}: deleted {deleted} BLAST DB file(s)')
            asm_filename = genome_json.get('assembly_fasta_file')
            if asm_filename:
                idx_path = os.path.join(genome.path, asm_filename + '.idx')
                if os.path.exists(idx_path):
                    os.remove(idx_path)
                    print(f'{genome_id}: deleted stale assembly FASTA index')
            continue

        if pre_check.has_pending_v3_files:
            names = ', '.join(os.path.basename(p) for p in pre_check.pending_files)
            print(f'{genome_id}: WARNING: partial upgrade detected ({names}). '
                  f'Use --promote to finish, or remove .v3 files manually and re-run to restart.')
            skipped_pending += 1
            continue

        gbk_filename = genome_json.get('cds_tool_gbk_file')
        if not gbk_filename:
            print(f'{genome_id}: SKIP: no cds_tool_gbk_file in genome.json')
            skipped_errors += 1
            continue
        gbk_path = os.path.join(genome.path, gbk_filename)
        if not os.path.exists(gbk_path):
            print(f'{genome_id}: SKIP: GBK not found: {gbk_path}')
            skipped_errors += 1
            continue

        if _locus_tags_need_rename(gbk_path, genome_id):
            detected = GenBankFile(gbk_path).detect_locus_tag_prefix()
            print(f'{genome_id}: SKIP: locus tags use prefix {detected!r}, expected {genome_id!r}. '
                  f'Fix manually with rename_genbank / rename_gff / rename_fasta / rename_eggnog '
                  f'before re-running the upgrade.')
            skipped_needs_rename += 1
            continue

        gbk_stem = os.path.splitext(gbk_path)[0]
        v3_created = set()   # every .v3 path touched: for cleanup on exception
        v3_to_orig = {}      # {v3_path: original_path}: for promotion

        try:
            # 2a. Rewrite GBK contig IDs → gbk.v3
            gbk_v3 = gbk_path + '.v3'
            v3_created.add(gbk_v3)
            contig_map = GenBankFile(gbk_path).normalize_contigs(
                out=gbk_v3, genome_id=genome_id, contig_format=contig_format)
            v3_to_orig[gbk_v3] = gbk_path
            print(f'{genome_id}: created {os.path.basename(gbk_v3)} ({len(contig_map)} contig IDs updated)')

            # 2b. Update assembly FNA contig headers
            if contig_map:
                asm_filename = genome_json.get('assembly_fasta_file')
                if asm_filename:
                    asm_path = os.path.join(genome.path, asm_filename)
                    if os.path.exists(asm_path):
                        asm_v3 = asm_path + '.v3'
                        v3_created.add(asm_v3)
                        renamed = _apply_contig_map_to_fna(asm_path, asm_v3, contig_map)
                        v3_to_orig[asm_v3] = asm_path
                        print(f'{genome_id}: created {os.path.basename(asm_v3)} ({renamed} contig headers updated)')
                        if renamed == 0:
                            print(f'{genome_id}: WARNING: assembly FNA: no contig headers matched '
                                  f'(seqid format may not be supported; plain or gnl|X|id expected). '
                                  f'Check {os.path.basename(asm_path)} manually.')

            # 2c. Rewrite GFF contig seqids.
            if not create_from_file and contig_map:
                gff_filename = genome_json.get('cds_tool_gff_file')
                gff_path = os.path.join(genome.path, gff_filename) if gff_filename else gbk_stem + '.gff'
                if os.path.exists(gff_path):
                    gff_v3 = gff_path + '.v3'
                    v3_created.add(gff_v3)
                    seqid_changed = _apply_contig_map_to_gff(gff_path, gff_v3, contig_map)
                    print(f'{genome_id}: created {os.path.basename(gff_v3)} ({seqid_changed} seqids updated)')
                    if seqid_changed == 0:
                        print(f'{genome_id}: WARNING: {os.path.basename(gff_path)}: no seqids matched '
                              f'(contig ID format in GFF may not be supported). Check manually.')
                    v3_to_orig[gff_v3] = gff_path


        except Exception as e:
            for v3_path in v3_created:
                try:
                    os.unlink(v3_path)
                except OSError:
                    pass
            if not isinstance(e, ValueError):
                logging.error('Unexpected error processing %s (%s)', genome_id, genome.path)
            raise

        # create_only: stop here; .v3 files are on disk, originals untouched.
        if create_only:
            names = ', '.join(os.path.basename(p) for p in sorted(v3_to_orig))
            print(f'{genome_id}: .v3 files created ({names}). Re-run with --promote to promote.')
            succeeded += 1
            continue

        # 3. Archive originals → tar.gz, move .v3 → originals.
        extra_backup = [gbk_stem + ext for ext in ('.fna', '.gff', '.faa', '.ffn')] if create_from_file else None
        archive = _promote_v3_files(v3_to_orig, genome_dir=genome.path, genome_id=genome_id,
                                    extra_backup=extra_backup)
        print(f'{genome_id}: archived originals → {os.path.basename(archive)}')

        # 3b. Regenerate gff/faa/ffn from the updated GBK (create_from_file only;
        #     otherwise the promoted .v3 files already contain the updated contig seqids).
        if create_from_file:
            gbk_final = GenBankFile(gbk_path)
            for ext, create_fn in [
                ('.fna', gbk_final.create_fna),
                ('.gff', gbk_final.create_gff),
                ('.faa', gbk_final.create_faa),
                ('.ffn', gbk_final.create_ffn),
            ]:
                out = gbk_stem + ext
                if os.path.exists(out):
                    os.remove(out)
                create_fn(out)

        # 4. Post-check
        post_check = check_genome_v3(genome.path, genome_id, deep=False, contig_format=contig_format)
        if not post_check.is_v3:
            reasons = '; '.join(post_check.issues)
            print(f'\n{"=" * 60}', file=sys.stderr)
            print(f'ERROR: {genome_id}: post-check failed', file=sys.stderr)
            for issue in post_check.issues:
                print(f'  - {issue}', file=sys.stderr)
            print(f'{"=" * 60}\n', file=sys.stderr)
            raise ValueError(f'{genome_id}: post-check failed: {reasons}')
        print(f'{genome_id}: done (post-check OK)')
        succeeded += 1

        # 5. Delete BLAST databases and stale sequence index files
        deleted = 0
        for fname in os.listdir(genome.path):
            if any(fname.endswith(ext) for ext in _BLAST_EXTENSIONS):
                os.remove(os.path.join(genome.path, fname))
                deleted += 1
        if deleted:
            print(f'{genome_id}: deleted {deleted} BLAST DB file(s)')
        asm_filename = genome_json.get('assembly_fasta_file')
        if asm_filename:
            idx_path = os.path.join(genome.path, asm_filename + '.idx')
            if os.path.exists(idx_path):
                os.remove(idx_path)
                print(f'{genome_id}: deleted stale assembly FASTA index')

    # ── Summary ──────────────────────────────────────────────────────────────
    if create_only:
        print(f'\nSummary: {succeeded} .v3 files created')
        print('Re-run with --promote to archive originals and promote .v3 files.')
        return

    def _extra(counts):
        parts = []
        if counts.get('rename'):
            parts.append(f'{counts["rename"]} skipped (wrong locus tag prefix; fix manually and re-run)')
        if counts.get('errors'):
            parts.append(f'{counts["errors"]} skipped (config/file errors)')
        return (', ' + ', '.join(parts)) if parts else ''

    blockers = skipped_needs_rename + skipped_errors

    if promote:
        print(f'\nSummary: {succeeded} promoted, {skipped_not_ready} not ready'
              + _extra({'rename': skipped_needs_rename, 'errors': skipped_errors}))
        if skipped_not_ready > 0:
            print(f'Not bumping to version 3: {skipped_not_ready} genome(s) still need --create_only. '
                  f'Run --create_only on remaining genomes, then re-run --promote.')
        elif blockers > 0:
            print(f'Not bumping to version 3: {blockers} genome(s) still need attention (see above).')
        else:
            set_folder_structure_version(new_version=3, folder_structure_dir=folder_structure_dir)
    else:
        print(f'\nSummary: {succeeded} migrated, {skipped_pending} skipped (pending .v3)'
              + _extra({'rename': skipped_needs_rename, 'errors': skipped_errors}))
        if skipped_pending > 0:
            print(f'Not bumping to version 3: {skipped_pending} genome(s) have pending .v3 files. '
                  f'Run --promote to finish them.')
        elif blockers > 0:
            print(f'Not bumping to version 3: {blockers} genome(s) still need attention (see above).')
        else:
            set_folder_structure_version(new_version=3, folder_structure_dir=folder_structure_dir)




def _gbk_to_assembly(gbk_path: str, fna_path: str) -> tuple[int, str]:
    """Regenerate assembly FNA from GBK sequences (overwrites fna_path). Returns (contig count, backup path)."""
    from Bio import SeqIO
    orig_descriptions = []
    with open(fna_path) as f:
        for line in f:
            if line.startswith('>'):
                parts = line[1:].strip().split(None, 1)
                orig_descriptions.append(parts[1] if len(parts) > 1 else '')
    backup = fna_path + '.bak'
    shutil.copy2(fna_path, backup)
    records = list(SeqIO.parse(gbk_path, 'genbank'))
    with open(fna_path, 'w') as f:
        for i, rec in enumerate(records):
            desc = orig_descriptions[i] if i < len(orig_descriptions) else ''
            f.write(f'>{rec.id} {desc}\n' if desc else f'>{rec.id}\n')
            seq = str(rec.seq)
            for j in range(0, len(seq), 60):
                f.write(seq[j:j + 60] + '\n')
    return len(records), backup


def check_assembly_compatibility(folder_structure_dir: str = None, skip_ignored: bool = False) -> None:
    """
    Check whether GBK and assembly FNA contig IDs match for each genome.
    Interactively offers to regenerate the assembly FNA from GBK sequences if they don't match.
    """
    folder_structure_dir = _get_folder_structure_dir(folder_structure_dir)
    ok = fixed = incompatible = skipped = 0

    for genome in loop_genomes(folder_structure_dir=folder_structure_dir, skip_ignored=skip_ignored):
        genome_id = genome.identifier
        genome_json = genome.json

        gbk_filename = genome_json.get('cds_tool_gbk_file')
        if not gbk_filename:
            print(f'{genome_id}: SKIP: no cds_tool_gbk_file in genome.json')
            skipped += 1
            continue
        gbk_path = os.path.join(genome.path, gbk_filename)
        if not os.path.exists(gbk_path):
            print(f'{genome_id}: SKIP: GBK not found: {gbk_path}')
            skipped += 1
            continue

        asm_filename = genome_json.get('assembly_fasta_file')
        if not asm_filename:
            print(f'{genome_id}: SKIP: no assembly_fasta_file in genome.json')
            skipped += 1
            continue
        fna_path = os.path.join(genome.path, asm_filename)
        if not os.path.exists(fna_path):
            print(f'{genome_id}: SKIP: assembly FNA not found: {fna_path}')
            skipped += 1
            continue

        gbk_ids = GenBankFile(gbk_path).get_contig_ids()
        gbk_id_set = set(gbk_ids)

        fna_ids = FastaFile(fna_path).get_contig_ids()
        fna_id_set = set(fna_ids)

        if fna_id_set == gbk_id_set:
            print(f'{genome_id}: OK')
            ok += 1
        elif ({_strip_gnl_prefix(i) for i in fna_id_set} == gbk_id_set
              or {_strip_version_suffix(i) for i in gbk_id_set} == fna_id_set
              or gbk_id_set == {_strip_version_suffix(i) for i in fna_id_set}):
            print(f'{genome_id}: OK')
            ok += 1
        else:
            while set(FastaFile(fna_path).get_contig_ids()) != gbk_id_set:
                fna_ids = FastaFile(fna_path).get_contig_ids()
                n_preview = 5
                print(f'{genome_id}: GBK and assembly contig IDs do not match '
                      f'(GBK: {len(gbk_ids)}, FNA: {len(fna_ids)}).')
                print(f'  FNA: {fna_path}')
                print(f'  {"GBK":30s}  {"FNA (current)":30s}')
                print(f'  {"---":30s}  {"---":30s}')
                for g, a in zip(gbk_ids[:n_preview], fna_ids[:n_preview]):
                    print(f'  {g:30s}  {a}')
                if max(len(gbk_ids), len(fna_ids)) > n_preview:
                    print(f'  ... ({max(len(gbk_ids), len(fna_ids)) - n_preview} more)')
                if len(gbk_ids) != len(fna_ids):
                    print(f'{genome_id}: contig counts differ; cannot auto-fix.')
                    response = input('  Fix manually, then (r)etry, or (s)kip? [r/s] ').strip().lower()
                    if response == 'r':
                        continue
                    break
                response = input('  (y) regenerate assembly from GBK, (r) retry after manual fix, (s) skip? [y/r/s] ').strip().lower()
                if response == 'y':
                    n, backup = _gbk_to_assembly(gbk_path, fna_path)
                    print(f'{genome_id}: assembly regenerated from GBK ({n} contigs); '
                          f'original backed up to {os.path.basename(backup)}')
                elif response == 'r':
                    continue
                else:
                    print(f'{genome_id}: skipped.')
                    break
            if set(FastaFile(fna_path).get_contig_ids()) == gbk_id_set:
                fixed += 1
            else:
                incompatible += 1

    print(f'\nSummary: {ok} OK, {fixed} fixed, {incompatible} incompatible, {skipped} skipped')


def check_v3(folder_structure_dir: str = None, genome_dir: str = None, genome_id: str = None,
             deep: bool = False, contig_format: str = '_scf{n}'):
    """
    Check v3 compatibility.

    Pass --folder_structure_dir to check all genomes, or --genome_dir for one genome.
    genome_id defaults to the basename of genome_dir.
    Use --deep to also check custom annotation files (default: GBK + assembly FNA only).
    """
    if genome_dir:
        if genome_id is None:
            genome_id = os.path.basename(genome_dir.rstrip('/'))
        result = check_genome_v3(genome_dir, genome_id, deep=deep, contig_format=contig_format)
        print(result.summary(genome_id))
        return

    folder_structure_dir = _get_folder_structure_dir(folder_structure_dir)
    ok = issues = 0
    for genome in loop_genomes(folder_structure_dir=folder_structure_dir):
        result = check_genome_v3(genome.path, genome.identifier, deep=deep, contig_format=contig_format)
        print(result.summary(genome.identifier))
        if result.has_pending_v3_files or not result.is_v3:
            issues += 1
        else:
            ok += 1
    print(f'\nSummary: {ok} OK, {issues} with issues')


def main():
    from fire import Fire

    Fire({
        'get_current_version': get_folder_structure_version,
        '1_to_2': from_1_to_2,
        '2_to_3': from_2_to_3,
        'check_v3': check_v3,
        'check_assembly_compatibility': check_assembly_compatibility,
    })


if __name__ == '__main__':
    main()
