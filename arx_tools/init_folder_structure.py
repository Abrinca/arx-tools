import os
import json
import shutil
from urllib import request
from .utils import PACKAGE_ROOT
from . import __folder_structure_version__


def download_go_data(out: str) -> None:
    source_url = 'https://current.geneontology.org/ontology/go.obo'

    print(f'Converting {source_url} -> {out}')

    def go_generator(io) -> [str]:
        go_entry = []

        line = io.readline()
        while line:
            if line == b'[Term]\n':
                yield go_entry
                go_entry.clear()

            go_entry.append(line.decode('utf-8'))
            line = io.readline()

        yield go_entry

    def get_name(entry: list) -> str:
        for line in entry:
            if line.startswith('name: '):
                return line.rstrip()[6:]
        raise TypeError(F'The go.obo file seems to have a wrong format! broken entry: {entry}')

    def get_go(entry: list) -> str:
        entry = entry[1]
        assert entry.startswith('id: GO:') and len(entry) == 15, f'Bad entry in go.obo: {entry}, len={len(entry)}'
        assert entry[7:14].isnumeric()
        return entry[4:14]

    with request.urlopen(source_url) as source_handle, open(out, 'w') as target_handle:
        gos = go_generator(io=source_handle)

        # skip first entry
        file_head = next(gos)
        assert not file_head[0].startswith(
            '[Term]'), F'The go.obo file seems to have a wrong format! file_head looks wrong: {file_head}'

        # save regular entries to file
        for entry in gos:
            target_handle.write(F'{get_go(entry)}\t{get_name(entry)}\n')


def download_kegg_data(src: str, out: str, remove_prefix: str = '', add_prefix: str = '') -> None:
    source_url = f'http://rest.kegg.jp/list/{src}'

    print(f'Converting {source_url} -> {out}')

    with request.urlopen(source_url) as source_handle, open(out, 'w') as target_handle:
        for line in source_handle:
            target_handle.write(f'{add_prefix}{line.decode("utf-8").removeprefix(remove_prefix)}')


def download_sl_data(out: str) -> None:
    # https://www.uniprot.org/locations -> Share -> Generate URL for API
    source_url = 'https://rest.uniprot.org/locations/stream?compressed=false' \
                 '&fields=id,name,definition&format=tsv&query=*'

    print(f'Converting {source_url} -> {out}')

    error_msg = 'UniProt must have changed its format. Please contact the developer. error={error}'

    try:
        source_handle = request.urlopen(source_url)
    except Exception:
        raise AssertionError(f'Failed to download {source_url}.\n{error_msg}')

    with open(out, 'w') as target_handle:
        first_line = source_handle.readline().decode('utf-8')
        assert first_line == 'Subcellular location ID\tName\tDescription\n', error_msg.format(error=first_line)

        for line in source_handle:
            line = line.decode('utf-8').strip().split('\t')
            assert len(line) == 3, error_msg.format(error=f'{len(line)=}; {line=}')
            sl, name, description = line

            target_handle.write(f'{sl}\t{name} ({description})\n')

    source_handle.close()

def init_folder_structure(folder_structure_dir: str = None) -> None:
    """
    Creates a basic OpenGenomeBrowser folders structure.

    Result:
        folder_structure
        ├── organisms
        ├── annotations.json
        ├── annotation-descriptions
        │   ├── SL.tsv
        │   ├── KO.tsv
        │   ├── KR.tsv
        │   ├── EC.tsv
        │   └── GO.tsv
        ├── orthologs
        └── pathway-maps
            ├── type_dictionary.json
            └── svg


    :param folder_structure_dir: Path to the root of the OpenGenomeBrowser folder structure. (Will contain 'organisms' folder.)
    """
    if folder_structure_dir is None:
        assert 'FOLDER_STRUCTURE' in os.environ, \
            f'Cannot find the folder_structure. ' \
            f'Please set --folder_structure_dir or environment variable FOLDER_STRUCTURE'
        folder_structure_dir = os.environ['FOLDER_STRUCTURE']

    assert os.path.isdir(
        os.path.dirname(folder_structure_dir)), f'Parent dir of {folder_structure_dir=} does not exist!'

    # make main dir
    os.makedirs(folder_structure_dir, exist_ok=True)

    # set version
    version_file = f'{folder_structure_dir}/version.json'
    if not os.path.exists(version_file):
        with open(version_file, 'w') as f:
            json.dump({'folder_structure_version': __folder_structure_version__}, f, indent=4)

    # make organisms dir (empty)
    os.makedirs(f'{folder_structure_dir}/organisms', exist_ok=True)

    # make orthologs dir (empty)
    os.makedirs(f'{folder_structure_dir}/orthologs', exist_ok=True)

    # make pathway maps dir and content
    os.makedirs(f'{folder_structure_dir}/pathway-maps/svg', exist_ok=True)
    type_dict = f'{folder_structure_dir}/pathway-maps/type_dictionary.json'
    if not os.path.exists(type_dict):
        with open(type_dict, 'w') as f:
            f.write('{}')

    # Create annotations.json
    annotations_dst = f'{folder_structure_dir}/annotations.json'
    if not os.path.exists(annotations_dst):
        shutil.copy(src=f'{PACKAGE_ROOT}/data/annotations.json', dst=annotations_dst)

    # download annotation descriptions (skip files that already exist)
    annotation_descriptions_dir = f'{folder_structure_dir}/annotation-descriptions'
    os.makedirs(annotation_descriptions_dir, exist_ok=True)
    if not os.path.exists(f'{annotation_descriptions_dir}/SL.tsv'):
        download_sl_data(out=f'{annotation_descriptions_dir}/SL.tsv')
    if not os.path.exists(f'{annotation_descriptions_dir}/KR.tsv'):
        download_kegg_data(src='rn', out=f'{annotation_descriptions_dir}/KR.tsv', remove_prefix='rn:')
    if not os.path.exists(f'{annotation_descriptions_dir}/KG.tsv'):
        download_kegg_data(src='ko', out=f'{annotation_descriptions_dir}/KG.tsv', remove_prefix='ko:')
    if not os.path.exists(f'{annotation_descriptions_dir}/EC.tsv'):
        download_kegg_data(src='enzyme', out=f'{annotation_descriptions_dir}/EC.tsv', remove_prefix='ec:', add_prefix='EC:')
    if not os.path.exists(f'{annotation_descriptions_dir}/GO.tsv'):
        download_go_data(out=f'{annotation_descriptions_dir}/GO.tsv')


def main():
    import fire

    fire.Fire(init_folder_structure)


if __name__ == '__main__':
    main()
