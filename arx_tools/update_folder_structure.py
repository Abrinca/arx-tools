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
        assert 'FOLDER_STRUCTURE' in os.environ, f'Cannot find the folder_structure. Please set --folder_structure_dir or environment variable FOLDER_STRUCTURE'
        folder_structure_dir = os.environ['FOLDER_STRUCTURE']
    assert os.path.isdir(folder_structure_dir), f'Could not find the folder_structure. Folder does not exist: {folder_structure_dir}'
    return folder_structure_dir


def set_folder_structure_version(new_version: int, folder_structure_dir: str) -> None:
    assert type(folder_structure_dir) is str
    version_file = f'{folder_structure_dir}/version.json'

    with open(version_file) as f:
        version_dict = json.load(f)

    version_dict['folder_structure_version'] = new_version

    with open(version_file, 'w') as f:
        json.dump(version_dict, f, indent=4)

    print()
    print(f'Successfully updated to folder structure version {new_version}!')


def ask(v_from: int, v_to: int, actions: [str], folder_structure_dir: str):
    assert type(folder_structure_dir) is str

    current_version = get_folder_structure_version(folder_structure_dir)
    assert current_version == v_from, \
        f'Cannot proceed: Folder structure version mismatch.\n' \
        f'This script expects version {v_from}, but folder_structure/version.json says version {current_version}.'

    question = f'Upgrade folder structure from version {v_from} to {v_to}:'
    for action in actions:
        question += f'\n - {action}'
    question += '\n\nProceed?'
    if not query_yes_no(question=question, default='yes'):
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
            print(f'{genome.identifier}: SKIPPED — permission denied on genome.json')
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

        COG = {}  # default

        eggnog_files = [f for f in genome_json['custom_annotations'] if f['type'].startswith('eggnog')]
        for file in eggnog_files:
            path = os.path.join(genome.path, file['file'])
            try:
                COG = EggnogFile(file=path).cog_categories()
            except AssertionError as e:
                logging.info(msg=str(e))
                pass

        print(f'{genome.identifier}: adding COG={COG}')
        genome_json['COG'] = COG
        genome.replace_json(genome_json)

    set_folder_structure_version(new_version=v_to, folder_structure_dir=folder_structure_dir)


def _apply_lt_map_to_file(src: str, dst: str, lt_map: dict, is_eggnog: bool = False) -> None:
    """
    Write a locus-tag-renamed copy of a tab-separated annotation file to dst.

    First column of each non-comment line is the locus_tag.
    Eggnog files may use "prefix|locus_tag" format — only the part after | is mapped.
    """
    with open(src) as f:
        lines = f.readlines()

    result = []
    for line in lines:
        if line.startswith('#') or not line.strip():
            result.append(line)
            continue
        cols = line.split('\t', 1)
        raw_tag = cols[0]
        if is_eggnog and '|' in raw_tag:
            prefix, locus_tag = raw_tag.rsplit('|', 1)
            new_tag = f'{prefix}|{lt_map.get(locus_tag, locus_tag)}'
        else:
            new_tag = lt_map.get(raw_tag, raw_tag)
        result.append(new_tag + ('\t' + cols[1] if len(cols) > 1 else '\n'))

    with open(dst, 'w') as f:
        f.writelines(result)


def _apply_lt_map_to_fasta(src: str, dst: str, lt_map: dict) -> int:
    """Write a locus-tag-renamed copy of a FASTA file to dst. Returns rename count."""
    with open(src) as f:
        lines = f.readlines()
    renamed = 0
    result = []
    for line in lines:
        if line.startswith('>'):
            parts = line[1:].split(None, 1)
            old_lt = parts[0]
            if old_lt in lt_map:
                rest = (' ' + parts[1]) if len(parts) > 1 else '\n'
                line = f'>{lt_map[old_lt]}{rest}'
                renamed += 1
        result.append(line)
    with open(dst, 'w') as f:
        f.writelines(result)
    return renamed


def _apply_contig_map_to_fna(src: str, dst: str, contig_map: dict) -> int:
    """Write a contig-header-renamed copy of a FASTA file to dst. Returns rename count."""
    with open(src) as f:
        lines = f.readlines()

    renamed = 0
    result = []
    for line in lines:
        if line.startswith('>'):
            parts = line[1:].split(None, 1)
            old_id = parts[0]
            if old_id in contig_map:
                rest = (' ' + parts[1]) if len(parts) > 1 else '\n'
                line = f'>{contig_map[old_id]}{rest}'
                renamed += 1
        result.append(line)

    with open(dst, 'w') as f:
        f.writelines(result)
    return renamed


_BLAST_EXTENSIONS = {
    # nucleotide v4/v5
    '.nhr', '.nin', '.nsq', '.nsi', '.nsd', '.ndb', '.not', '.ntf', '.nto', '.njs',
    # protein v4/v5
    '.phr', '.pin', '.psq', '.psi', '.psd', '.pdb', '.pot', '.ptf', '.pto', '.pjs',
    # alias files
    '.nal', '.pal',
}


def _promote_v3_files(v3_to_orig: dict, genome_dir: str, genome_id: str) -> str:
    """
    Back up original files into a single {genome_id}_v2_backup.tar.gz, then promote .v3 → originals.

    v3_to_orig maps {v3_path: original_path}. Files are archived with paths relative to genome_dir.
    Returns the path of the backup archive.
    """
    archive_path = os.path.join(genome_dir, f'{genome_id}_v2_backup.tar.gz')
    with tarfile.open(archive_path, 'w:gz') as tar:
        for orig_path in v3_to_orig.values():
            if os.path.exists(orig_path):
                tar.add(orig_path, arcname=os.path.relpath(orig_path, genome_dir))

    for v3_path, orig_path in v3_to_orig.items():
        if os.path.exists(orig_path):
            os.unlink(orig_path)
        shutil.move(v3_path, orig_path)

    return archive_path


def from_2_to_3(folder_structure_dir: str = None, skip_ignored=False, contig_format: str = '_scf{n}'):
    """
    Upgrade folder structure from v2 to v3.

    Per genome:
      1. Shallow v3 check — skip if already v3; warn if a partial upgrade (.v3 files) exists.
      2. Generate .v3 intermediate files for source files (gbk, assembly fna, annotations).
         Originals are untouched during this phase.
      3. Only once all .v3 files are successfully written: archive originals into a single
         {genome_id}_v2_backup.tar.gz and rename .v3 files into place.
      3b. Regenerate derived files (.fna, .gff, .faa, .ffn) from the promoted GBK.
      4. Post-check to verify success.
      5. Delete BLAST databases (they reference stale contig/locus IDs).
    """
    folder_structure_dir = _get_folder_structure_dir(folder_structure_dir)

    warnings.filterwarnings('ignore', message='.*malformed locus line.*', module='Bio')

    ask(
        v_from=2, v_to=3,
        actions=[
            'shallow-check each genome; skip if already v3',
            'generate .v3 intermediate files for source files (gbk, assembly fna, annotations)',
            'archive originals into {genome_id}_v2_backup.tar.gz and promote .v3 files into place',
            'post-check each genome to verify',
            'delete BLAST databases (will be rebuilt on next import)',
        ],
        folder_structure_dir=folder_structure_dir,
    )

    for genome in loop_genomes(folder_structure_dir=folder_structure_dir, skip_ignored=skip_ignored):
        genome_json = genome.json
        genome_id = genome.identifier

        # 1. Shallow check
        pre = check_genome_v3(genome.path, genome_id, deep=False, contig_format=contig_format)
        if pre.is_v3:
            print(f'{genome_id}: already v3, skipping')
            continue
        if pre.has_pending_v3_files:
            names = ', '.join(os.path.basename(p) for p in pre.pending_files)
            print(f'{genome_id}: WARNING — partial upgrade detected ({names}). '
                  f'Remove .v3 files manually and re-run to retry.')
            continue

        gbk_filename = genome_json.get('cds_tool_gbk_file')
        if not gbk_filename:
            print(f'{genome_id}: no cds_tool_gbk_file in genome.json, skipping')
            continue

        gbk_path = os.path.join(genome.path, gbk_filename)
        if not os.path.exists(gbk_path):
            print(f'{genome_id}: GBK not found at {gbk_path}, skipping')
            continue

        base = os.path.splitext(gbk_path)[0]
        v3_created = set()   # every .v3 path touched — for cleanup on failure
        v3_to_orig = {}      # {v3_path: original_path} — only successful files, for promotion
        failed = False

        # 2a. Normalize GBK → gbk.v3
        gbk_v3 = gbk_path + '.v3'
        v3_created.add(gbk_v3)
        try:
            contig_map, lt_map = GenBankFile(gbk_path).normalize(
                out=gbk_v3, genome_id=genome_id, contig_format=contig_format)
            v3_to_orig[gbk_v3] = gbk_path
            print(f'{genome_id}: created {os.path.basename(gbk_v3)} ({len(lt_map)} locus tags renamed)')
        except Exception as e:
            print(f'{genome_id}: ERROR normalizing GBK: {e}')
            failed = True

        # (derived files .fna/.gff/.faa/.ffn are regenerated after promotion — no .v3 needed)

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
                    except Exception as e:
                        print(f'{genome_id}: ERROR updating assembly FNA: {e}')
                        failed = True

        # 2d. Update custom annotation files
        if not failed and lt_map:
            ca_count = 0
            for ca in genome_json.get('custom_annotations', []):
                ca_path = os.path.join(genome.path, ca['file'])
                if not os.path.exists(ca_path):
                    print(f'{genome_id}: custom annotation not found: {ca["file"]}, skipping')
                    continue
                ca_v3 = ca_path + '.v3'
                v3_created.add(ca_v3)
                try:
                    _apply_lt_map_to_file(ca_path, ca_v3, lt_map, is_eggnog=ca['type'].startswith('eggnog'))
                    v3_to_orig[ca_v3] = ca_path
                    ca_count += 1
                except Exception as e:
                    print(f'{genome_id}: ERROR updating {ca["file"]}: {e}')
                    failed = True
            if ca_count:
                print(f'{genome_id}: created .v3 for {ca_count} annotation file(s)')

        # On failure: clean up every .v3 file we may have created, leave originals untouched
        if failed:
            for v3_path in v3_created:
                try:
                    os.unlink(v3_path)
                except OSError:
                    pass
            print(f'{genome_id}: FAILED — original files untouched. Fix errors and re-run.')
            continue

        # 3. Archive originals → tar.gz, move .v3 → originals
        archive = _promote_v3_files(v3_to_orig, genome_dir=genome.path, genome_id=genome_id)
        print(f'{genome_id}: archived originals → {os.path.basename(archive)}')

        # 3b. Update derived files: rename IDs/locus-tags in-place to preserve headers;
        #     regenerate .gff (fully structured, no free-form content).
        gbk_final = GenBankFile(gbk_path)
        for ext, apply_fn, create_fn in [
            ('.fna', lambda p: _apply_contig_map_to_fna(p, p + '.new', contig_map), gbk_final.create_fna),
            ('.faa', lambda p: _apply_lt_map_to_fasta(p, p + '.new', lt_map), gbk_final.create_faa),
            ('.ffn', lambda p: _apply_lt_map_to_fasta(p, p + '.new', lt_map), gbk_final.create_ffn),
        ]:
            out = base + ext
            try:
                if os.path.exists(out):
                    apply_fn(out)
                    os.replace(out + '.new', out)
                else:
                    create_fn(out)
            except Exception as e:
                print(f'{genome_id}: ERROR updating {ext}: {e}')
        gff_path = base + '.gff'
        if os.path.exists(gff_path):
            os.remove(gff_path)
        try:
            gbk_final.create_gff(gff_path)
        except Exception as e:
            print(f'{genome_id}: ERROR regenerating .gff: {e}')

        # 4. Post-check
        post = check_genome_v3(genome.path, genome_id, deep=False, contig_format=contig_format)
        if post.is_v3:
            print(f'{genome_id}: done (post-check OK)')
        else:
            print(f'{genome_id}: WARNING — post-check failed: {"; ".join(post.issues)}')

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
