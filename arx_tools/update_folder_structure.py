import os
import json
import logging
import shutil
import tempfile

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
        if genome.has_json:
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


def _apply_lt_map_to_annotation_file(path: str, lt_map: dict, is_eggnog: bool = False) -> None:
    """
    Rewrite a tab-separated annotation file using a full locus_tag rename map.

    First column of each non-comment line is the locus_tag.
    Eggnog files may use "prefix|locus_tag" format — only the part after | is mapped.
    Writes atomically via a temp file.
    """
    with open(path) as f:
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

    with tempfile.NamedTemporaryFile('w', delete=False, dir=os.path.dirname(path)) as tmp:
        tmp.writelines(result)
        tmp_path = tmp.name
    shutil.move(tmp_path, path)


_BLAST_EXTENSIONS = {
    # nucleotide v4/v5
    '.nhr', '.nin', '.nsq', '.nsi', '.nsd', '.ndb', '.not', '.ntf', '.nto', '.njs',
    # protein v4/v5
    '.phr', '.pin', '.psq', '.psi', '.psd', '.pdb', '.pot', '.ptf', '.pto', '.pjs',
    # alias files
    '.nal', '.pal',
}


def from_2_to_3(folder_structure_dir: str = None, skip_ignored=False):
    """
    Upgrade folder structure from v2 to v3.

    Per genome:
      - Backs up the .gbk file to {genome}.gbk.v2.bak
      - Canonicalizes contig IDs ({genome}_scf00001, …) and locus_tags ({genome}_00001, …)
      - Regenerates .fna, .gff, .faa, .ffn from the normalized .gbk
      - Deletes BLAST databases (they reference old contig/locus IDs)
    """
    folder_structure_dir = _get_folder_structure_dir(folder_structure_dir)

    ask(
        v_from=2, v_to=3,
        actions=[
            'back up each .gbk to {genome}.gbk.v2.bak',
            'rename contigs to {genome}_scf00001, … and locus_tags to {genome}_00001, …',
            'apply locus_tag rename map to all custom annotation files (eggnog, GO, EC, …)',
            'regenerate .fna, .gff, .faa, .ffn from the normalized .gbk',
            'delete BLAST databases (will be rebuilt on next import)',
        ],
        folder_structure_dir=folder_structure_dir,
    )

    for genome in loop_genomes(folder_structure_dir=folder_structure_dir, skip_ignored=skip_ignored):
        genome_json = genome.json
        gbk_filename = genome_json.get('cds_tool_gbk_file')
        if not gbk_filename:
            print(f'{genome.identifier}: no cds_tool_gbk_file in genome.json, skipping')
            continue

        gbk_path = os.path.join(genome.path, gbk_filename)
        if not os.path.exists(gbk_path):
            print(f'{genome.identifier}: .gbk not found at {gbk_path}, skipping')
            continue

        genome_id = genome.identifier
        base = os.path.splitext(gbk_path)[0]

        # 1. Backup
        backup_path = gbk_path + '.v2.bak'
        if not os.path.exists(backup_path):
            shutil.copy2(gbk_path, backup_path)
            print(f'{genome_id}: backed up .gbk → {os.path.basename(backup_path)}')

        # 2. Normalize .gbk (rename contigs + locus_tags) via temp file
        try:
            with tempfile.NamedTemporaryFile(delete=False, suffix='.gbk') as tmp:
                tmp_path = tmp.name
            _, lt_map = GenBankFile(gbk_path).normalize(out=tmp_path, genome_id=genome_id)
            shutil.move(tmp_path, gbk_path)
            print(f'{genome_id}: normalized .gbk ({len(lt_map)} locus tags renamed)')
        except Exception as e:
            print(f'{genome_id}: ERROR normalizing .gbk: {e}')
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
            continue

        # 3. Rename locus tags in custom annotation files
        custom_annotations = genome_json.get('custom_annotations', [])
        for ca in custom_annotations:
            ca_path = os.path.join(genome.path, ca['file'])
            if not os.path.exists(ca_path):
                print(f'{genome_id}: custom annotation file not found: {ca["file"]}, skipping')
                continue
            try:
                _apply_lt_map_to_annotation_file(ca_path, lt_map, is_eggnog=ca['type'].startswith('eggnog'))
            except Exception as e:
                print(f'{genome_id}: ERROR renaming locus tags in {ca["file"]}: {e}')
        if custom_annotations:
            print(f'{genome_id}: updated locus tags in {len(custom_annotations)} custom annotation file(s)')

        gbk = GenBankFile(gbk_path)

        # 5. Regenerate derived files
        for ext, create_fn in [
            ('.fna', gbk.create_fna),
            ('.gff', gbk.create_gff),
            ('.faa', gbk.create_faa),
            ('.ffn', gbk.create_ffn),
        ]:
            derived = base + ext
            try:
                if os.path.exists(derived):
                    os.remove(derived)
                create_fn(derived)
            except Exception as e:
                print(f'{genome_id}: ERROR regenerating {ext}: {e}')

        # 6. Delete BLAST databases
        deleted = 0
        for fname in os.listdir(genome.path):
            if any(fname.endswith(ext) for ext in _BLAST_EXTENSIONS):
                os.remove(os.path.join(genome.path, fname))
                deleted += 1
        if deleted:
            print(f'{genome_id}: deleted {deleted} BLAST DB file(s)')

        print(f'{genome_id}: done')

    set_folder_structure_version(new_version=3, folder_structure_dir=folder_structure_dir)


def main():
    from fire import Fire

    Fire({
        'get_current_version': get_folder_structure_version,
        '1_to_2': from_1_to_2,
        '2_to_3': from_2_to_3,
    })


if __name__ == '__main__':
    main()
