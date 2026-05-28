import os
import json
import logging
import shutil
import tarfile
import warnings

from .check_v3 import check_genome_v3
from .folder_looper import FolderLooper, FolderGenome
from .rename_eggnog import EggnogFile
from .rename_genbank import GenBankFile
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


def _apply_gene_tag_map_to_file(src: str, dst: str, gene_tag_map: dict, is_eggnog: bool = False) -> tuple[int, int]:
    """
    Write a locus-tag-renamed copy of a tab-separated annotation file to dst.

    First column of each non-comment line is the locus_tag.
    Eggnog files may use "prefix|locus_tag" format: only the part after | is mapped.

    Returns (found, total): counts of data lines where the locus tag was / was not in gene_tag_map.
    """
    with open(src) as f:
        lines = f.readlines()

    matched = total = 0
    result = []
    for line in lines:
        if line.startswith('#') or not line.strip():
            result.append(line)
            continue
        total += 1
        cols = line.split('\t', 1)
        raw_tag = cols[0]
        if is_eggnog and '|' in raw_tag:
            prefix, locus_tag = raw_tag.rsplit('|', 1)
            if locus_tag in gene_tag_map:
                matched += 1
                new_tag = f'{prefix}|{gene_tag_map[locus_tag]}'
            else:
                new_tag = raw_tag
        else:
            if raw_tag in gene_tag_map:
                matched += 1
                new_tag = gene_tag_map[raw_tag]
            else:
                new_tag = raw_tag
        result.append(new_tag + ('\t' + cols[1] if len(cols) > 1 else '\n'))

    with open(dst, 'w') as f:
        f.writelines(result)
    return matched, total


def _extend_gene_tag_map(gene_tag_map: dict) -> dict:
    """
    Return a copy of gene_tag_map extended with 5-digit-padded variant keys.

    Annotation files and derived FASTA files (FAA/FFN) may have been generated
    against arx-assigned 5-digit locus tags (e.g. GENOME_ID_00001) even when the
    source GBK stores external locus tags such as NCBI RefSeq IDs.  In those cases
    gene_tag_map only maps NCBI_TAG → GENOME_ID_000001 and a direct lookup of
    GENOME_ID_00001 fails.  This function adds GENOME_ID_00001 → GENOME_ID_000001
    entries derived from the v3 values already in gene_tag_map so those files can
    be renamed correctly.
    """
    extended = dict(gene_tag_map)
    for v3_tag in gene_tag_map.values():
        sep = v3_tag.rfind('_')
        if sep < 0:
            continue
        digits_str = v3_tag[sep + 1:]
        if not digits_str.isdigit() or len(digits_str) < 6:
            continue
        number = int(digits_str)
        if number > 99999:
            continue
        five_digit_key = v3_tag[:sep + 1] + f'{number:05d}'
        if five_digit_key not in extended:
            extended[five_digit_key] = v3_tag
    return extended


def _apply_gene_tag_map_to_fasta(src: str, dst: str, gene_tag_map: dict) -> tuple[int, int]:
    """
    Write a locus-tag-renamed copy of a FASTA file to dst.

    Returns (renamed, total): count of headers renamed and total header count.
    Callers can compare the two to warn when no renames occurred despite a non-trivial map.
    """
    with open(src) as f:
        lines = f.readlines()
    renamed = total = 0
    result = []
    for line in lines:
        if line.startswith('>'):
            total += 1
            parts = line[1:].split(None, 1)
            old_gene_tag = parts[0]
            if old_gene_tag in gene_tag_map:
                rest = (' ' + parts[1]) if len(parts) > 1 else '\n'
                line = f'>{gene_tag_map[old_gene_tag]}{rest}'
                renamed += 1
        result.append(line)
    with open(dst, 'w') as f:
        f.writelines(result)
    return renamed, total


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
                raise ValueError(f"No contig id found in {line!r}")
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


def _rename_gff_feature_line(line: str, contig_map: dict, gene_tag_map: dict) -> tuple[str, int, int]:
    """Rename seqid and attributes in a GFF3 feature line; return (line, seqid_changed, attr_renamed)."""
    cols = line.split('\t')
    if len(cols) != 9:
        return line, 0, 0

    seqid_changed = 0
    new_contig = _resolve_contig_id(cols[0], contig_map)
    if new_contig is not None:
        cols[0] = new_contig
        seqid_changed = 1

    trailing = '\n' if cols[8].endswith('\n') else ''
    attr_renamed = 0
    new_attrs = []
    for attr in cols[8].rstrip('\n').split(';'):
        if '=' in attr:
            key, val = attr.split('=', 1)
            new_val = _lookup_gene_tag(val, gene_tag_map)
            new_attrs.append(f'{key}={new_val}')
            # Only ID=/Parent= renames require format-specific logic; locus_tag= always succeeds directly.
            if key in ('ID', 'Parent') and new_val != val:
                attr_renamed += 1
        else:
            new_attrs.append(attr)
    cols[8] = ';'.join(new_attrs) + trailing
    return '\t'.join(cols), seqid_changed, attr_renamed


def _apply_maps_to_gff(src: str, dst: str, contig_map: dict, gene_tag_map: dict) -> tuple[int, int]:
    """
    Write a renamed copy of a GFF file to dst.

    Returns (seqid_changed, attr_renamed). A zero attr_renamed despite a non-empty
    gene_tag_map means the GFF uses an unsupported feature-ID format.
    """
    with open(src) as f:
        lines = f.readlines()

    seqid_changed = 0
    attr_renamed = 0
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
            line, sc, ar = _rename_gff_feature_line(line, contig_map, gene_tag_map)
            seqid_changed += sc
            attr_renamed += ar
        result.append(line)

    with open(dst, 'w') as f:
        f.writelines(result)
    return seqid_changed, attr_renamed


def _resolve_contig_id(contig_id: str, contig_map: dict) -> str | None:
    """Look up contig_id, stripping Prokka's gnl|X| prefix if needed (GBK keys are bare IDs)."""
    new_id = contig_map.get(contig_id)
    if new_id is not None:
        return new_id
    if '|' in contig_id:
        base = contig_id.rsplit('|', 1)[1]
        if base in contig_map:
            return contig_map[base]
    return None


# Prokka appends feature-type suffixes to locus tags in GFF3 ID attributes (e.g. GENOME_00001_gene).
_PROKKA_FEATURE_SUFFIXES = (
    '_gene', '_tRNA', '_rRNA', '_tmRNA', '_repeat_region',
    '_misc_RNA', '_CRISPR', '_ncRNA',
)

# PGAP prefixes feature-type to locus tags in GFF3 IDs (e.g. gene-RS00001, rna-RS00150).
# cds- may use a protein accession (e.g. cds-WP_012211193.1); those aren't in gene_tag_map and pass through.
_PGAP_FEATURE_PREFIXES = ('gene-', 'cds-', 'rna-', 'exon-', 'id-')


def _lookup_gene_tag(val: str, gene_tag_map: dict) -> str:
    """Look up val with fallbacks for Prokka suffix and PGAP prefix conventions; return val unchanged if not found."""
    new_val = gene_tag_map.get(val)
    if new_val is not None:
        return new_val
    # Prokka: locus_tag_feature_type (trailing underscore-prefixed suffix)
    for suffix in _PROKKA_FEATURE_SUFFIXES:
        if val.endswith(suffix):
            base = val[:-len(suffix)]
            if base in gene_tag_map:
                return gene_tag_map[base] + suffix
    # PGAP: feature_type-locus_tag (leading hyphen-terminated prefix)
    for prefix in _PGAP_FEATURE_PREFIXES:
        if val.startswith(prefix):
            rest = val[len(prefix):]
            if rest in gene_tag_map:
                return prefix + gene_tag_map[rest]
            # exon-locus_tag-N: locus tag followed by a hyphen + integer exon number
            if '-' in rest:
                base, sep, exon_n = rest.rpartition('-')
                if exon_n.isdigit() and base in gene_tag_map:
                    return f'{prefix}{gene_tag_map[base]}{sep}{exon_n}'
    return val


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

    Per genome:
      1. Shallow v3 check: skip if already v3; warn if a partial upgrade (.v3 files) exists.
      2. Generate .v3 intermediate files for source files (gbk, assembly fna, annotations).
         Originals are untouched during this phase.
      3. Only once all .v3 files are successfully written: archive originals into a single
         {genome_id}_v2_backup.tar.gz and rename .v3 files into place.
      3b. Update derived files (.fna, .gff, .faa, .ffn) by renaming IDs in-place to preserve
          custom content. With --create_from_file, regenerate all derived files from the GBK.
      4. Post-check to verify success.
      5. Delete BLAST databases (they reference stale contig/locus IDs).

    Two-step workflow:
      --create_only  Stop after step 2 — generate .v3 files without promoting them.
                     Inspect the generated files, then re-run with --promote to finish.
      --promote      Skip step 2 — promote .v3 files left by a previous --create_only run.
    """
    folder_structure_dir = _get_folder_structure_dir(folder_structure_dir)

    warnings.filterwarnings('ignore', message='.*malformed locus line.*', module='Bio')
    warnings.filterwarnings('ignore', message='.*Premature end of file.*', module='Bio')
    warnings.filterwarnings('ignore', message='.*Expected sequence length.*', module='Bio')

    derived_action = (
        'regenerate derived files (.fna, .gff, .faa, .ffn) from GBK (--create_from_file)'
        if create_from_file else
        'generate .v3 for existing derived files (.fna, .gff, .faa, .ffn); promote all .v3 files together'
    )
    if create_only:
        actions = [
            'shallow-check each genome; skip if already v3',
            'generate .v3 intermediate files for source files (gbk, assembly fna, annotations)',
            derived_action.split(';')[0].strip() if not create_from_file else derived_action,
            '(promotion skipped — re-run with --promote to archive originals and promote)',
        ]
    elif promote:
        actions = [
            'shallow-check each genome; skip if already v3 or no pending .v3 files',
            'archive originals into {genome_id}_v2_backup.tar.gz and promote pending .v3 files into place',
            'post-check each genome to verify',
            'delete BLAST databases (will be rebuilt on next import)',
        ]
    else:
        actions = [
            'shallow-check each genome; skip if already v3',
            'generate .v3 intermediate files for source files (gbk, assembly fna, annotations)',
            'archive originals into {genome_id}_v2_backup.tar.gz and promote .v3 files into place',
            derived_action,
            'post-check each genome to verify',
            'delete BLAST databases (will be rebuilt on next import)',
        ]
    ask(v_from=2, v_to=3, actions=actions, folder_structure_dir=folder_structure_dir)
    genomes_iter = loop_genomes(folder_structure_dir=folder_structure_dir, skip_ignored=skip_ignored)

    succeeded = 0
    failed_count = 0
    post_check_failed = []
    skipped_not_ready = 0   # --promote: genomes still v2 with no pending .v3 files
    skipped_pending = 0     # normal run: genomes with leftover .v3 files (warn and skip)

    for genome in loop_genomes(folder_structure_dir=folder_structure_dir, skip_ignored=skip_ignored):
        genome_json = genome.json
        genome_id = genome.identifier

        # 1. Shallow check
        pre_check = check_genome_v3(genome.path, genome_id, deep=False, contig_format=contig_format)
        if pre_check.is_v3:
            print(f'{genome_id}: already v3, skipping')
            continue

        if promote:
            # --promote: skip creation, only handle .v3 files left by a previous --create_only run.
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
            if post_check.is_v3:
                print(f'{genome_id}: done (post-check OK)')
                succeeded += 1
            else:
                reasons = '; '.join(post_check.issues)
                print(f'{genome_id}: WARNING: post-check failed: {reasons}')
                post_check_failed.append(genome_id)
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
            print(f'{genome_id}: no cds_tool_gbk_file in genome.json, skipping')
            continue

        gbk_path = os.path.join(genome.path, gbk_filename)
        if not os.path.exists(gbk_path):
            print(f'{genome_id}: GBK not found at {gbk_path}, skipping')
            continue

        gbk_stem = os.path.splitext(gbk_path)[0]
        v3_created = set()   # every .v3 path touched: for cleanup on failure
        v3_to_orig = {}      # {v3_path: original_path}: only successful files, for promotion
        failed = False

        # 2a. Normalize GBK → gbk.v3
        gbk_v3 = gbk_path + '.v3'
        v3_created.add(gbk_v3)
        try:
            contig_map, gene_tag_map = GenBankFile(gbk_path).normalize(
                out=gbk_v3, genome_id=genome_id, contig_format=contig_format)
            v3_to_orig[gbk_v3] = gbk_path
            print(f'{genome_id}: created {os.path.basename(gbk_v3)} ({len(gene_tag_map)} locus tags renamed)')
        except Exception as e:
            print(f'{genome_id}: ERROR normalizing GBK: {e}')
            failed = True

        # 2b. Build extended map so annotation/derived files that use arx-assigned
        #     5-digit locus tags (e.g. GENOME_ID_00001) are also covered when the
        #     GBK stores external/NCBI locus tags.
        if not failed:
            extended_gene_tag_map = _extend_gene_tag_map(gene_tag_map)
            # True when at least one locus tag actually changes (not an identity map).
            # Used below to suppress false-positive "no renames" warnings for genomes
            # whose locus tags are already v3-compliant.
            _nontrivial_lt_map = any(k != v for k, v in gene_tag_map.items())

        # 2c. Update assembly FNA contig headers
        if not failed and contig_map:
            asm_filename = genome_json.get('assembly_fasta_file')
            if asm_filename:
                asm_path = os.path.join(genome.path, asm_filename)
                if os.path.exists(asm_path):
                    asm_v3 = asm_path + '.v3'
                    v3_created.add(asm_v3)
                    try:
                        renamed = _apply_contig_map_to_fna(asm_path, asm_v3, contig_map)
                        v3_to_orig[asm_v3] = asm_path
                        print(f'{genome_id}: created {os.path.basename(asm_v3)} ({renamed} contig headers updated)')
                        if renamed == 0:
                            print(f'{genome_id}: WARNING: assembly FNA: no contig headers matched '
                                  f'— the seqid format may not be supported (plain or gnl|X|id). '
                                  f'Check {os.path.basename(asm_path)} manually.')
                    except Exception as e:
                        print(f'{genome_id}: ERROR updating assembly FNA: {e}')
                        failed = True

        # 2d. Update custom annotation files
        # Continue loop on failure to report all errors before giving up.
        if not failed and gene_tag_map:
            annotation_count = 0
            for annotation in genome_json.get('custom_annotations', []):
                annotation_path = os.path.join(genome.path, annotation['file'])
                if not os.path.realpath(annotation_path).startswith(os.path.realpath(genome.path) + os.sep):
                    print(f'{genome_id}: WARNING: annotation path escapes genome dir, skipping: {annotation["file"]}')
                    continue
                if not os.path.exists(annotation_path):
                    print(f'{genome_id}: custom annotation not found: {annotation["file"]}, skipping')
                    continue
                ann_type = annotation['type']
                is_eggnog = ann_type.startswith('eggnog')
                annotation_v3 = annotation_path + '.v3'
                v3_created.add(annotation_v3)
                try:
                    found, total = _apply_gene_tag_map_to_file(annotation_path, annotation_v3, extended_gene_tag_map, is_eggnog=is_eggnog)
                    if total > 0 and found == 0:
                        print(f'{genome_id}: WARNING: {annotation["file"]} (type={ann_type!r}): '
                              f'no locus tags matched: file may be in an unsupported format or use a different column')
                    v3_to_orig[annotation_v3] = annotation_path
                    annotation_count += 1
                except Exception as e:
                    print(f'{genome_id}: ERROR updating {annotation["file"]}: {e}')
                    failed = True
            if annotation_count:
                print(f'{genome_id}: created .v3 for {annotation_count} annotation file(s)')

        # 2e. Create .v3 for existing derived files (skipped when create_from_file, as they'll be regenerated).
        if not failed and not create_from_file:
            for ext in ('.fna', '.gff', '.faa', '.ffn'):
                orig_path = gbk_stem + ext
                if not os.path.exists(orig_path):
                    continue
                v3_path = orig_path + '.v3'
                v3_created.add(v3_path)
                try:
                    if ext == '.fna':
                        renamed = _apply_contig_map_to_fna(orig_path, v3_path, contig_map)
                        print(f'{genome_id}: created {os.path.basename(v3_path)} ({renamed} contig headers updated)')
                        if contig_map and renamed == 0:
                            print(f'{genome_id}: WARNING: {os.path.basename(orig_path)}: no contig headers '
                                  f'matched — seqid format may not be supported. Check manually.')
                    elif ext == '.gff':
                        seqid_changed, attr_renamed = _apply_maps_to_gff(
                            orig_path, v3_path, contig_map, extended_gene_tag_map)
                        print(f'{genome_id}: created {os.path.basename(v3_path)} '
                              f'({seqid_changed} seqids, {attr_renamed} attributes updated)')
                        if contig_map and seqid_changed == 0:
                            print(f'{genome_id}: WARNING: {os.path.basename(orig_path)}: no seqids '
                                  f'matched — contig ID format in GFF may not be supported. Check manually.')
                        if _nontrivial_lt_map and attr_renamed == 0:
                            print(f'{genome_id}: WARNING: {os.path.basename(orig_path)}: no attribute '
                                  f'values renamed — ID=/Parent= format may not be supported '
                                  f'(known: Prokka _gene/_tRNA suffixes, PGAP gene-/cds-/rna-/exon- prefixes). '
                                  f'Check manually.')
                    else:  # .faa / .ffn
                        renamed, total = _apply_gene_tag_map_to_fasta(
                            orig_path, v3_path, extended_gene_tag_map)
                        print(f'{genome_id}: created {os.path.basename(v3_path)} ({renamed}/{total} headers updated)')
                        if _nontrivial_lt_map and total > 0 and renamed == 0:
                            print(f'{genome_id}: WARNING: {os.path.basename(orig_path)}: no FASTA '
                                  f'headers matched — locus tag format may differ from what is expected. '
                                  f'Check manually.')
                    v3_to_orig[v3_path] = orig_path
                except Exception as e:
                    print(f'{genome_id}: ERROR creating {os.path.basename(v3_path)}: {e}')
                    failed = True

        # On failure: clean up every .v3 file we may have created, leave originals untouched
        if failed:
            for v3_path in v3_created:
                try:
                    os.unlink(v3_path)
                except OSError:
                    pass
            print(f'{genome_id}: FAILED: original files untouched. Fix errors and re-run.')
            failed_count += 1
            continue

        # create_only: stop here — .v3 files are on disk, originals untouched.
        if create_only:
            names = ', '.join(os.path.basename(p) for p in sorted(v3_to_orig))
            print(f'{genome_id}: .v3 files created ({names}). Re-run with --promote to promote.')
            succeeded += 1
            continue

        # 3. Archive originals → tar.gz, move .v3 → originals.
        #    For create_from_file=True, also back up derived files as extras (no .v3 intermediate, will be regenerated).
        extra_backup = [gbk_stem + ext for ext in ('.fna', '.gff', '.faa', '.ffn')] if create_from_file else None
        archive = _promote_v3_files(v3_to_orig, genome_dir=genome.path, genome_id=genome_id,
                                    extra_backup=extra_backup)
        print(f'{genome_id}: archived originals → {os.path.basename(archive)}')

        # 3b. For create_from_file=True: regenerate all derived files from the normalized GBK.
        #     For create_from_file=False: derived .v3 files were already promoted in step 3.
        create_from_file_error = False
        if create_from_file:
            gbk_final = GenBankFile(gbk_path)
            for ext, create_fn in [
                ('.fna', gbk_final.create_fna),
                ('.gff', gbk_final.create_gff),
                ('.faa', gbk_final.create_faa),
                ('.ffn', gbk_final.create_ffn),
            ]:
                out = gbk_stem + ext
                try:
                    if os.path.exists(out):
                        os.remove(out)
                    create_fn(out)
                except Exception as e:
                    print(f'{genome_id}: ERROR creating {ext}: {e}')
                    create_from_file_error = True

        # 4. Post-check
        post_check = check_genome_v3(genome.path, genome_id, deep=False, contig_format=contig_format)
        if post_check.is_v3 and not create_from_file_error:
            print(f'{genome_id}: done (post-check OK)')
            succeeded += 1
        else:
            reasons = '; '.join(post_check.issues)
            if create_from_file_error:
                reasons = ('derived file regeneration error(s)' + ('; ' + reasons if reasons else '')).rstrip('; ')
            print(f'{genome_id}: WARNING: post-check failed: {reasons}')
            post_check_failed.append(genome_id)

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
        print(f'\nSummary: {succeeded} .v3 files created, {failed_count} failed')
        if failed_count == 0:
            print('Re-run with --promote to archive originals and promote .v3 files.')
        return

    if promote:
        print(f'\nSummary: {succeeded} promoted, {skipped_not_ready} not ready, '
              f'{failed_count} failed, {len(post_check_failed)} failed post-check')
    else:
        print(f'\nSummary: {succeeded} migrated, {skipped_pending} skipped (pending .v3), '
              f'{failed_count} failed, {len(post_check_failed)} failed post-check')

    if failed_count > 0 or post_check_failed:
        if post_check_failed:
            print(f'Not bumping to version 3: post-check failures: {", ".join(post_check_failed)}. Fix and re-run.')
        else:
            print('Not bumping to version 3: fix errors and re-run.')
    elif skipped_not_ready > 0:
        print(f'Not bumping to version 3: {skipped_not_ready} genome(s) still need --create_only. '
              f'Run --create_only on remaining genomes, then re-run --promote.')
    elif skipped_pending > 0:
        print(f'Not bumping to version 3: {skipped_pending} genome(s) have pending .v3 files. '
              f'Run --promote to finish them.')
    else:
        set_folder_structure_version(new_version=3, folder_structure_dir=folder_structure_dir)


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
    })


if __name__ == '__main__':
    main()
