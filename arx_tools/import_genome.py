import logging
import os
import json
import yaml
import shutil
import tempfile
from glob import glob
from textwrap import shorten
from typing import Union
from schema import SchemaError

from . import __folder_structure_version__
from .utils import entrez_organism_to_taxid, GenomeFile, merge_json, get_folder_structure_version, WorkingDirectory
from .rename_genbank import GenBankFile
from .rename_gff import GffFile
from .rename_fasta import FastaFile
from .rename_eggnog import EggnogFile
from .parse_busco import parse_busco
from .rename_custom_annotations import CustomAnnotationFile
from .metadata_schemas import \
    organism_json_schema, genome_json_dummy, genome_json_schema, organism_json_dummy


class ImportException(Exception):
    pass


class ImportSettings:
    settings: dict[str:str]
    default_settings = {
        'organism_template': {},
        'genome_template': {},
        'import_actions': [
            {'type': 'copy', 'from': '*', 'to': '{original_path}', 'expected': True},
        ],
        'file_finder': {
            'fna': {'glob': '*.fna', 'expected': False},
            'gbk': {'glob': '*.gbk', 'expected': 1},
            'gff': {'glob': '*.gff', 'expected': False},
            'faa': {'glob': '*.faa', 'expected': False},
            'sqn': {'glob': '*.sqn', 'expected': False},
            'ffn': {'glob': '*.ffn', 'expected': False},
            'eggnog': {'glob': '*.emapper.annotations', 'expected': False},
            'yaml': {'glob': '*.yaml', 'expected': False},
            'busco': {'glob': '*_busco.txt', 'expected': False},
            'custom_annotations': [
                {'glob': f'*.{anno_type}', 'anno_type': anno_type, 'expected': False}
                for anno_type in ('GC', 'GP', 'EP', 'ED', 'EO', 'EC', 'KG', 'KR', 'GO', 'SL', 'OL')]
        }
    }

    def __init__(self, settings: Union[dict, str] = None):
        if settings is None:
            settings = os.environ.get('ARX_IMPORT_SETTINGS', None)

        if type(settings) is dict:
            pass
        elif type(settings) is str:
            with open(settings) as f:
                settings = json.load(f)
        else:
            settings = {}

        self.settings = self.default_settings | settings  # PEP-584: dict union: overwrite defaults with new settings
        self.settings['file_finder'] = self.default_settings['file_finder'] | settings.get('file_finder', {})

        assert set(self.settings.keys()) == set(self.default_settings.keys()), \
            f'ARX_IMPORT_SETTINGS must contain these JSON keys: {set(self.default_settings.keys())}! ' \
            f'reality: {self.settings.keys()}'

    @staticmethod
    def _format_path(path: str, genome: str, organism: str, src: str = None) -> str:
        if src is None:
            return path.format(
                genome=genome,
                organism=organism,
                assembly=genome.rsplit('.', 1)[0]
            )
        else:
            basename = os.path.basename(src)
            return path.format(
                original_path=src,
                basename=basename,
                suffix=basename.rsplit('.', 1)[-1] if '.' in basename else '',
                genome=genome,
                organism=organism,
                assembly=genome.rsplit('.', 1)[0]
            )

    @staticmethod
    def _copy(src: str, dst: str):
        if os.path.exists(dst):
            logging.warning(f'Overwriting: {src} -> {dst}')
        copy_fn = shutil.copy2 if os.path.isfile(src) else shutil.copytree
        os.makedirs(os.path.dirname(dst), exist_ok=True)  # create parent dir if nonexistent
        copy_fn(src=src, dst=dst)

    @classmethod
    def copy(cls, source_dir: str, target_dir: str, genome: str, organism: str, action: dict):
        from_ = action['from']
        to = action['to']
        expected = action.get('expected', True)

        with WorkingDirectory(source_dir):
            files = glob(from_)
            cls.check_expected(files, expected, from_)
            for src in files:
                rel_dst = cls._format_path(to, genome, organism, src)
                dst = os.path.join(target_dir, rel_dst)
                if os.path.isdir(dst):
                    logging.warning(f'Overwriting directory: {src} >>{action}>> {rel_dst}')
                    shutil.rmtree(dst)
                elif os.path.isfile(dst):
                    logging.warning(f'Overwriting file: {src} >>{action}>> {rel_dst}')
                    os.remove(dst)
                else:
                    logging.info(f'{src} >>{action}>> {rel_dst}')
                cls._copy(src=src, dst=dst)

    @classmethod
    def link(cls, target_dir: str, genome: str, organism: str, action: dict):
        from_ = cls._format_path(action['from'], genome, organism)
        to = cls._format_path(action['to'], genome, organism)
        expected = action.get('expected', True)
        assert type(expected) is bool, f'Failed to execute link action: "expected" must be true or false! {action=}'

        dst = os.path.join(target_dir, to)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        os.symlink(src=from_, dst=dst)
        if expected:
            assert os.path.exists(dst), f'Failed to execute link action: destination {dst=} does not exist! {action=}'

    def execute_actions(self, source_dir: str, target_dir: str, genome: str, organism: str) -> None:
        for action in self.settings['import_actions']:
            action_type = action['type']
            if action_type == 'copy':
                self.copy(source_dir, target_dir, genome, organism, action)
            elif action_type == 'link':
                self.link(target_dir, genome, organism, action)
            else:
                raise AssertionError(f'Could not execute action: type must be "copy". {action=}')

    @staticmethod
    def check_expected(files: [str], expected: Union[None, bool, int], glob_pattern: str):
        if type(expected) is int:
            if len(files) != expected:
                raise ImportException(f'Error: {files=} glob={glob_pattern}\n'
                                      f'Found {len(files)}, not {expected} files!')
        elif expected is True:
            if len(files) == 0:
                raise ImportException(f'Error: Found no files using glob={glob_pattern}!')
        else:
            if expected is not False:
                raise ImportException(f'Error: config is bad. type must be integer or boolean. '
                                      f'{expected=} {type(expected)=}')

    def find_files(self, type_: str, root_dir: str) -> [str]:
        settings = self.settings['file_finder'][type_]
        glob_pattern = settings['glob']
        expected = settings.get('expected', False)

        with WorkingDirectory(root_dir):
            files = glob(glob_pattern)

        logging.info(f'Found {len(files)} files of type={type_} using glob={glob_pattern}')
        self.check_expected(files, expected, glob_pattern)
        return files

    def find_file(self, type_: str, root_dir: str, as_class=None, expected: bool = True) -> Union[str, GenomeFile, None]:
        files = self.find_files(type_, root_dir)

        if len(files) == 1:
            abs_path = os.path.join(root_dir, files[0])
            if as_class is None:
                return abs_path
            else:
                return as_class(abs_path)
        else:
            if expected:
                raise AssertionError(
                    f'Error: found {len(files)} files of {type_=}: {files=}')
            else:
                logging.info(f'Found no {type_} files.')
                return None

    def find_custom_annotations(self, root_dir: str):
        annotations = []

        eggnog_file = self.find_file(type_='eggnog', root_dir=root_dir, as_class=EggnogFile, expected=False)
        if eggnog_file:
            annotations.append(eggnog_file)

        with WorkingDirectory(root_dir):
            for custom_annotation in self.settings['file_finder']['custom_annotations']:
                glob_pattern = custom_annotation['glob']
                files = glob(glob_pattern)
                expected = custom_annotation.get('expected', False)
                assert len(files) < 2, f'Found multiple {custom_annotation["anno_type"]}: {files=}'

                if files:
                    annotations.append(CustomAnnotationFile(
                        file=os.path.join(root_dir, files[0]),
                        custom_annotation_type=custom_annotation['anno_type'])
                    )
                if not files and expected:
                    raise ImportException(f'Error: Found no custom-file using glob={glob_pattern}!')
        return annotations


def autodetect_organism_genome(root_dir: str) -> (str, str):
    with WorkingDirectory(root_dir):
        gbks = glob('*.gbk')
        last_error = None
        for gbk in gbks:
            try:
                strain, locus_tag_prefix = GenBankFile(file=gbk).detect_strain_locus_tag_prefix()
                organism, genome = strain, locus_tag_prefix.rstrip('_')
                logging.info(f'autodetected from gbk: {organism=} {genome=}')
                return organism, genome
            except Exception as e:
                last_error = e
    cause = f': {last_error}' if last_error else ''
    raise AssertionError(
        f'Failed to automatically detect organism and genome name from {gbks}{cause}. '
        f'Please specify the names manually.'
    )


def rename_all(root_dir: str, gbk: GenBankFile, files: [GenomeFile], new_prefix: str, old_prefix: str = None):
    if not old_prefix:
        old_prefix = gbk.detect_locus_tag_prefix()

    assert new_prefix != old_prefix, \
        f'old and new locus_tag_prefix are the same! {old_prefix=} {new_prefix=}'

    kwargs = dict(new_locus_tag_prefix=new_prefix, old_locus_tag_prefix=old_prefix, update_path=False)

    with tempfile.TemporaryDirectory() as rename_tempdir, WorkingDirectory(root_dir):
        for file in files:
            temp_file = os.path.join(rename_tempdir, 'tempfile')
            file.rename(out=temp_file, **kwargs)
            os.remove(file.path)
            os.rename(src=temp_file, dst=file.path)


def load_yaml_metadata(submol_yaml: str) -> (dict, dict):
    organism_yaml, genome_yaml = {}, {}
    with open(submol_yaml) as f:
        submol_yaml = yaml.safe_load(f)

    if 'organism' in submol_yaml and 'genus_species' in submol_yaml['organism']:
        organism_yaml['taxid'] = entrez_organism_to_taxid(submol_yaml['organism']['genus_species'])

    if 'biosample' in submol_yaml:
        genome_yaml['biosample_accession'] = submol_yaml['biosample']
    if 'bioproject' in submol_yaml:
        genome_yaml['bioproject_accession'] = submol_yaml['bioproject']

    if 'publications' in submol_yaml and len(submol_yaml['publications']):
        genome_yaml['literature_references'] = [{
            'url': f"https://pubmed.ncbi.nlm.nih.gov/{p['publication']['pmid']}/",
            'name': shorten(p['publication']['title'], width=30)}
            for p in submol_yaml['publications']
        ]

    return organism_yaml, genome_yaml


def load_cog_metadata(custom_annotations: [GenomeFile]) -> dict:
    for file in custom_annotations:
        if type(file) is EggnogFile:
            try:
                cog = file.cog_categories()
                return {'COG': cog}
            except AssertionError as e:
                logging.info(f'Failed to extract COG information from {file.path}. {str(e)}')
                pass
    return {}  # not eggnog file


def add_files_to_json(genome_json: dict, files: dict, custom_annotations) -> dict:
    def relname(key):
        file = files[key]
        return None if file is None else os.path.basename(file.path)

    genome_json['cds_tool_faa_file'] = relname('faa')
    genome_json['cds_tool_ffn_file'] = relname('ffn')
    genome_json['cds_tool_gbk_file'] = relname('gbk')
    genome_json['cds_tool_gff_file'] = relname('gff')
    genome_json['cds_tool_sqn_file'] = relname('sqn')
    genome_json['assembly_fasta_file'] = relname('fna')
    genome_json['custom_annotations'] = [
        {'date': ca.date_str(), 'file': os.path.basename(ca.path), 'type': ca.custom_annotation_type}
        for ca in custom_annotations
    ]
    return genome_json

def gather_metadata(import_settings: ImportSettings, root_dir: str, files: [GenomeFile],
                    custom_annotations: [GenomeFile], organism_dir: str, import_dir: str,
                    organism: str, genome: str):
    '''
    Load metadata from:
      - pgap_submol.yaml
      - *.gbk
      - *_busco.txt
      - organism.json and genome.json
    :return:
    '''

    # start with dummy jsons
    organism_json = organism_json_dummy.copy()
    genome_json = genome_json_dummy.copy()

    # add import_settings
    organism_json.update(import_settings.settings['organism_template'])
    genome_json.update(import_settings.settings['genome_template'])

    # add pgap_submol.yaml
    try:
        organism_yaml, genome_yaml = load_yaml_metadata(import_settings.find_file(type_='yaml', root_dir=root_dir))
        organism_json.update(organism_yaml)
        genome_json.update(genome_yaml)
    except AssertionError as e:
        logging.info(f'Failed to load metadata from yaml: {e}')

    # add *.gbk
    organism_gbk, genome_gbk = files['gbk'].metadata()
    organism_json.update(organism_gbk)
    genome_json.update(genome_gbk)

    # add _busco.txt
    try:
        busco_file = import_settings.find_file(type_='busco', root_dir=root_dir)
        genome_json['BUSCO'] = parse_busco(busco_file)
    except AssertionError:
        pass

    # add COG from eggnog
    genome_json.update(load_cog_metadata(custom_annotations))

    # add organism.json from folder structure
    organism_json = merge_json(organism_json, os.path.join(organism_dir, 'organism.json'))

    # add organism.json / genome.json from import_dir
    organism_json = merge_json(organism_json, os.path.join(import_dir, 'organism.json'))
    genome_json = merge_json(genome_json, os.path.join(import_dir, 'genome.json'))
    genome_json.pop('contig_format', None)  # input-only config, not stored in output

    # add elementary identifiers
    organism_json['name'] = organism
    organism_json['representative'] = genome
    genome_json['identifier'] = genome

    # add files
    genome_json = add_files_to_json(genome_json, files, custom_annotations)

    # validate metadata files
    try:
        organism_json_schema.validate(organism_json)
    except SchemaError as e:
        logging.warning(f'FAILED TO CREATE A VALID organism.json! {str(e)}')
        raise e

    try:
        genome_json_schema.validate(genome_json)
    except SchemaError as e:
        logging.warning(f'FAILED TO CREATE A VALID genome.json! {str(e)}')
        raise e

    return organism_json, genome_json


def check_files_(genome_id: str, locus_tag_prefix: str, files: dict, custom_annotations: [GenomeFile],
                 contig_format: str = '_scf{n}') -> None:
    files['gbk'].validate_contig_ids(genome_id=genome_id, contig_format=contig_format)
    files['gbk'].validate_locus_tags(locus_tag_prefix=locus_tag_prefix)
    if files['fna'] is not None:
        files['fna'].validate_contig_ids(genome_id=genome_id, contig_format=contig_format)
    files['gff'].validate_locus_tags(locus_tag_prefix=locus_tag_prefix)
    files['faa'].validate_locus_tags(locus_tag_prefix=locus_tag_prefix)
    files['ffn'].validate_locus_tags(locus_tag_prefix=locus_tag_prefix)
    for ca in custom_annotations:
        ca.validate_locus_tags(locus_tag_prefix=locus_tag_prefix)


def import_genome(
        import_dir: str,
        folder_structure_dir: str = None,
        organism: str = None,
        genome: str = None,
        rename: bool = False,
        check_files: bool = True,
        import_settings: str = None,
        pause: bool = False
):
    """
    Easily import files into OpenGenomeBrowser folder structure.

    :param import_dir: Folder with files to import. Required: [.fna, .faa, .gbk, .gff] Optional: [.ffn, .sqn, custom-annotation-files]
    :param folder_structure_dir: Path to the root of the OpenGenomeBrowser folder structure. (Must contain 'organisms' folder.)
    :param organism: Name of the organism.
    :param genome: Identifier of the genome. Must start with organism. May be identical to organism.
    :param rename: Locus tag prefixes must match the genome identifier. If this is not the case, this script can automatically rename relevant files.
    :param check_files: If true, check if locus tag prefixes match genome identifier.
    :param import_settings: Path to import settings file. Alternatively, set the environment variable ARX_IMPORT_SETTINGS.
    :param pause: Wait after import_actions / before file_finder
    """
    import_dir = os.path.abspath(import_dir)

    if folder_structure_dir is None:
        assert 'FOLDER_STRUCTURE' in os.environ, \
            f'Cannot find the folder_structure. ' \
            f'Please set --folder_structure_dir or environment variable FOLDER_STRUCTURE'
        folder_structure_dir = os.environ['FOLDER_STRUCTURE']

    folder_structure_dir = os.path.abspath(folder_structure_dir)

    organisms_dir = f'{folder_structure_dir}/organisms'
    assert os.path.isdir(organisms_dir), f'Cannot import files: {organisms_dir=} does not exist.'

    current_folder_structure_version = get_folder_structure_version(folder_structure_dir)
    assert current_folder_structure_version == __folder_structure_version__, \
        f'Before importing any genomes, the folder structure needs to be updated to match OpenGenomeBrowser Tools.\n' \
        f'Current version: {current_folder_structure_version}, expected: {__folder_structure_version__}\n' \
        f'Use the script update_folder_structure perform the upgrade!'

    assert os.path.isdir(import_dir), f'Cannot import files: {import_dir=} does not exist.'

    # Read optional contig_format from import_dir/genome.json before any processing.
    contig_format = '_scf{n}'
    _import_gj_path = os.path.join(import_dir, 'genome.json')
    if os.path.isfile(_import_gj_path):
        with open(_import_gj_path) as _f:
            _import_gj = json.load(_f)
        if 'contig_format' in _import_gj:
            contig_format = _import_gj['contig_format']

    import_settings = ImportSettings(import_settings)

    if organism is None or genome is None:
        _organism, _genome = autodetect_organism_genome(import_dir)
        if organism is None:
            organism = _organism
        if genome is None:
            genome = _genome

    # genome names can consist of integers -.-
    organism, genome = str(organism), str(genome)

    if not genome.startswith(organism):
        raise ImportException(
            f'Genome identifier ({genome!r}) must start with organism name ({organism!r}). '
            f'Consider using {organism!r}_{genome!r} instead.'
        )

    organism_dir = os.path.join(organisms_dir, organism)
    genome_dir = os.path.join(organism_dir, 'genomes', genome)
    assert not os.path.exists(genome_dir), f'Could not import {organism}:{genome}: {genome_dir=} already exists!'

    with tempfile.TemporaryDirectory() as work_dir:
        import_settings.execute_actions(import_dir, work_dir, genome, organism)

        if pause:
            print(f'Files are prepared here: {work_dir} Press enter to continue with import. Press Ctrl+C to abort.')
            input()

        gbk: GenBankFile = import_settings.find_file('gbk', root_dir=work_dir, as_class=GenBankFile)
        base = os.path.splitext(os.path.basename(gbk.path))[0]

        if rename:
            logging.info('Normalizing locus tags and contig IDs.')

            # If an fna is already present, match contigs by ID (not position) so the
            # canonical numbering follows GBK order and the fna gets renamed to match.
            _pre_fna = import_settings.find_file('fna', root_dir=work_dir, as_class=FastaFile, expected=False)
            if _pre_fna is not None:
                fna_contig_ids = _pre_fna.get_contig_ids()
                gbk_contig_ids = gbk.get_contig_ids()
                # Prokka writes FNA headers as "gnl|X|bare_id" (e.g. "gnl|C|ALNJDMAK_1")
                # while the GBK LOCUS line is just the bare id ("ALNJDMAK_1").
                # Strip the gnl|X| prefix before comparing so Prokka genomes pass.
                fna_contig_ids_norm = [cid.rsplit('|', 1)[1] if '|' in cid else cid for cid in fna_contig_ids]
                assert set(fna_contig_ids_norm) == set(gbk_contig_ids), (
                    f'FNA and GBK contain different contig IDs.\n'
                    f'  FNA only: {set(fna_contig_ids_norm) - set(gbk_contig_ids)}\n'
                    f'  GBK only: {set(gbk_contig_ids) - set(fna_contig_ids_norm)}'
                )
                # canonical IDs numbered by GBK order (consistent with the no-fna path)
                contig_id_map = {
                    gbk_id: f'{genome}{contig_format.format(n=i + 1)}'
                    for i, gbk_id in enumerate(gbk_contig_ids)
                }
                tmp_fna = _pre_fna.path + '.renaming'
                _pre_fna.rename_contig_ids(out=tmp_fna, new_ids=[contig_id_map[i] for i in fna_contig_ids_norm], update_path=False)
                os.replace(tmp_fna, _pre_fna.path)
                canonical_ids = [contig_id_map[i] for i in gbk_contig_ids]
            else:
                canonical_ids = None  # normalize generates them from contig_format

            tmp_path = gbk.path + '.normalizing'
            _, lt_map = gbk.normalize(out=tmp_path, genome_id=genome, contig_ids=canonical_ids, contig_format=contig_format)
            os.replace(tmp_path, gbk.path)

            # After normalization, any provided gff/faa/ffn have stale locus tags.
            # Remove them so they are regenerated from the normalized GBK below.
            with WorkingDirectory(work_dir):
                for _stale_pattern in ('*.gff', '*.faa', '*.ffn'):
                    for _stale_file in glob(_stale_pattern):
                        logging.info(f'Removing stale derived file after normalization: {_stale_file}')
                        os.remove(_stale_file)

            # Rename custom annotations (eggnog, .GC, etc.) using the exact locus tag map.
            for _ca in import_settings.find_custom_annotations(work_dir):
                tmp_ca = _ca.path + '.renaming'
                _ca.rename_by_map(out=tmp_ca, lt_map=lt_map, update_path=False)
                os.replace(tmp_ca, _ca.path)
                logging.info(f'Renamed locus tags in custom annotation: {_ca.path}')

        fna = import_settings.find_file('fna', root_dir=work_dir, as_class=FastaFile, expected=False)
        if fna is None:
            logging.info('Generating .fna from .gbk.')
            fna_path = os.path.join(work_dir, base + '.fna')
            gbk.create_fna(fna=fna_path)
            fna = FastaFile(fna_path)

        gff = import_settings.find_file('gff', root_dir=work_dir, as_class=GffFile, expected=False)
        if gff is None:
            logging.info('Generating .gff from .gbk.')
            gff_path = os.path.join(work_dir, base + '.gff')
            gbk.create_gff(gff=gff_path)
            gff = GffFile(gff_path)

        ffn = import_settings.find_file('ffn', root_dir=work_dir, as_class=FastaFile, expected=False)
        if ffn is None:
            logging.info('Generating .ffn from .gbk.')
            ffn_path = os.path.join(work_dir, base + '.ffn')
            gbk.create_ffn(ffn=ffn_path)
            ffn = FastaFile(ffn_path)

        faa = import_settings.find_file('faa', root_dir=work_dir, as_class=FastaFile, expected=False)
        if faa is None:
            logging.info('Generating .faa from .gbk.')
            faa_path = os.path.join(work_dir, base + '.faa')
            gbk.create_faa(faa=faa_path)
            faa = FastaFile(faa_path)

        sqn: GenomeFile = import_settings.find_file(
            'sqn', root_dir=work_dir, as_class=GenomeFile, expected=False)

        files = dict(fna=fna, gbk=gbk, ffn=ffn, faa=faa, gff=gff, sqn=sqn)

        custom_annotations = import_settings.find_custom_annotations(work_dir)

        organism_json, genome_json = gather_metadata(
            import_settings,
            root_dir=work_dir,
            files=files,
            custom_annotations=custom_annotations,
            organism_dir=organism_dir,
            import_dir=import_dir,
            organism=organism,
            genome=genome
        )

        if check_files:
            check_files_(genome_id=genome, locus_tag_prefix=f'{genome}_', files=files, custom_annotations=custom_annotations, contig_format=contig_format)

        # final movement
        os.makedirs(os.path.dirname(genome_dir), exist_ok=True)
        shutil.copytree(src=work_dir, dst=genome_dir, symlinks=True)

    with open(os.path.join(organism_dir, 'organism.json'), 'w') as f:
        json.dump(organism_json, f, indent=4)
    with open(os.path.join(genome_dir, 'genome.json'), 'w') as f:
        json.dump(genome_json, f, indent=4)


def main():
    import fire

    fire.Fire(import_genome)


if __name__ == '__main__':
    main()
