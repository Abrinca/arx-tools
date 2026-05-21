import json
import os
import tempfile
from unittest import TestCase

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

from arx_tools.check_v3 import check_genome_v3, V3CheckResult

GENOME_ID = 'FAM1079-i1-1.1'
CONTIG_FORMAT = '_scf{n}'


def _write_gbk(path: str, contigs: list[tuple[str, list[str]]]) -> None:
    """Write a minimal GenBank file. contigs: [(contig_id, [locus_tag, ...]), ...]"""
    records = []
    for contig_id, locus_tags in contigs:
        rec = SeqRecord(Seq('ATCGATCGATCG' * 5), id=contig_id, name=contig_id[:16], description='')
        rec.annotations['molecule_type'] = 'DNA'
        for i, lt in enumerate(locus_tags):
            feat = SeqFeature(FeatureLocation(i * 9, i * 9 + 9, strand=1), type='CDS')
            feat.qualifiers['locus_tag'] = [lt]
            rec.features.append(feat)
        records.append(rec)
    with open(path, 'w') as f:
        SeqIO.write(records, f, 'genbank')


def _write_fna(path: str, contigs: list[str]) -> None:
    """Write a minimal FASTA file with given contig IDs."""
    with open(path, 'w') as f:
        for contig_id in contigs:
            f.write(f'>{contig_id}\nATCGATCGATCG\n')


def _write_annotation(path: str, rows: list[tuple[str, str]]) -> None:
    """Write a minimal tab-separated annotation file."""
    with open(path, 'w') as f:
        f.write('# comment line\n')
        for locus_tag, value in rows:
            f.write(f'{locus_tag}\t{value}\n')


def _write_eggnog(path: str, rows: list[tuple[str, str, str]]) -> None:
    """Write a minimal eggnog-style file with prefix|locus_tag in first column."""
    with open(path, 'w') as f:
        f.write('# eggnog comment\n')
        for prefix, locus_tag, value in rows:
            f.write(f'{prefix}|{locus_tag}\t{value}\n')


def _setup_genome(tmp_dir: str, gbk_contigs=None, fna_contigs=None,
                  annotations=None, eggnog_rows=None) -> dict:
    """
    Create a minimal genome folder in tmp_dir and return the genome.json dict.

    gbk_contigs: [(contig_id, [locus_tag, ...]), ...]
    fna_contigs: [contig_id, ...]  (for assembly FNA)
    annotations: [(locus_tag, value), ...]  (for a KG annotation file)
    eggnog_rows: [(prefix, locus_tag, value), ...]
    """
    genome_json = {'identifier': GENOME_ID}

    if gbk_contigs is not None:
        gbk_path = os.path.join(tmp_dir, f'{GENOME_ID}.gbk')
        _write_gbk(gbk_path, gbk_contigs)
        genome_json['cds_tool_gbk_file'] = f'{GENOME_ID}.gbk'

    if fna_contigs is not None:
        fna_path = os.path.join(tmp_dir, f'{GENOME_ID}.fna')
        _write_fna(fna_path, fna_contigs)
        genome_json['assembly_fasta_file'] = f'{GENOME_ID}.fna'

    custom_annotations = []
    if annotations is not None:
        ca_path = os.path.join(tmp_dir, f'{GENOME_ID}.KG')
        _write_annotation(ca_path, annotations)
        custom_annotations.append({'file': f'{GENOME_ID}.KG', 'type': 'KG'})

    if eggnog_rows is not None:
        eg_path = os.path.join(tmp_dir, 'query_seqs.fa.emapper.annotations')
        _write_eggnog(eg_path, eggnog_rows)
        custom_annotations.append({'file': 'query_seqs.fa.emapper.annotations', 'type': 'eggnog'})

    if custom_annotations:
        genome_json['custom_annotations'] = custom_annotations

    with open(os.path.join(tmp_dir, 'genome.json'), 'w') as f:
        json.dump(genome_json, f)

    return genome_json


def _v3_gbk_contigs(n_contigs=2, n_loci=3):
    """Return v3-compliant GBK contigs for GENOME_ID."""
    return [
        (f'{GENOME_ID}_scf{i + 1}', [f'{GENOME_ID}_{str(j + 1).zfill(6)}' for j in range(n_loci)])
        for i in range(n_contigs)
    ]


def _v3_fna_contigs(n=2):
    return [f'{GENOME_ID}_scf{i + 1}' for i in range(n)]


class TestCheckV3Shallow(TestCase):
    def test_already_v3(self):
        with tempfile.TemporaryDirectory() as tmp:
            _setup_genome(tmp,
                          gbk_contigs=_v3_gbk_contigs(),
                          fna_contigs=_v3_fna_contigs())
            result = check_genome_v3(tmp, GENOME_ID)
            self.assertTrue(result.is_v3)
            self.assertFalse(result.has_pending_v3_files)
            self.assertEqual(result.issues, [])

    def test_bad_locus_tags(self):
        with tempfile.TemporaryDirectory() as tmp:
            _setup_genome(tmp,
                          gbk_contigs=[(f'{GENOME_ID}_scf1', ['OLD_00001', 'OLD_00002'])],
                          fna_contigs=_v3_fna_contigs(1))
            result = check_genome_v3(tmp, GENOME_ID)
            self.assertFalse(result.is_v3)
            self.assertTrue(any('locus_tag' in issue for issue in result.issues))

    def test_bad_contig_ids_in_gbk(self):
        with tempfile.TemporaryDirectory() as tmp:
            _setup_genome(tmp,
                          gbk_contigs=[('old_contig_1', [f'{GENOME_ID}_000001'])],
                          fna_contigs=[f'{GENOME_ID}_scf1'])
            result = check_genome_v3(tmp, GENOME_ID)
            self.assertFalse(result.is_v3)
            self.assertTrue(any('contig' in issue for issue in result.issues))

    def test_bad_assembly_fna_headers(self):
        with tempfile.TemporaryDirectory() as tmp:
            _setup_genome(tmp,
                          gbk_contigs=_v3_gbk_contigs(1),
                          fna_contigs=['old_contig_1'])  # non-v3 assembly FNA
            result = check_genome_v3(tmp, GENOME_ID)
            self.assertFalse(result.is_v3)
            self.assertTrue(any('FNA' in issue for issue in result.issues))

    def test_five_digit_locus_tags_not_v3(self):
        """zfill(5) output from the old normalizer must not pass as v3."""
        with tempfile.TemporaryDirectory() as tmp:
            _setup_genome(tmp,
                          gbk_contigs=[(f'{GENOME_ID}_scf1', [f'{GENOME_ID}_00001'])],
                          fna_contigs=[f'{GENOME_ID}_scf1'])
            result = check_genome_v3(tmp, GENOME_ID)
            self.assertFalse(result.is_v3)

    def test_seven_digit_locus_tags_accepted(self):
        """Genomes with >999,999 genes produce 7-digit tags via zfill(6): these must pass."""
        with tempfile.TemporaryDirectory() as tmp:
            _setup_genome(tmp,
                          gbk_contigs=[(f'{GENOME_ID}_scf1', [f'{GENOME_ID}_1000001'])],
                          fna_contigs=[f'{GENOME_ID}_scf1'])
            result = check_genome_v3(tmp, GENOME_ID)
            self.assertTrue(result.is_v3)

    def test_missing_genome_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = check_genome_v3(tmp, GENOME_ID)
            self.assertFalse(result.is_v3)
            self.assertIn('genome.json not found', result.issues)

    def test_missing_gbk_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            genome_json = {'identifier': GENOME_ID, 'cds_tool_gbk_file': 'missing.gbk'}
            with open(os.path.join(tmp, 'genome.json'), 'w') as f:
                json.dump(genome_json, f)
            result = check_genome_v3(tmp, GENOME_ID)
            self.assertFalse(result.is_v3)
            self.assertTrue(any('GBK not found' in issue for issue in result.issues))

    def test_no_gbk_in_json(self):
        """Genome without cds_tool_gbk_file: only FNA is checked."""
        with tempfile.TemporaryDirectory() as tmp:
            _setup_genome(tmp, fna_contigs=_v3_fna_contigs())
            result = check_genome_v3(tmp, GENOME_ID)
            self.assertTrue(result.is_v3)

    def test_summary_v3_ok(self):
        with tempfile.TemporaryDirectory() as tmp:
            _setup_genome(tmp, gbk_contigs=_v3_gbk_contigs(), fna_contigs=_v3_fna_contigs())
            result = check_genome_v3(tmp, GENOME_ID)
            self.assertIn('OK', result.summary(GENOME_ID))

    def test_summary_not_v3(self):
        with tempfile.TemporaryDirectory() as tmp:
            _setup_genome(tmp, gbk_contigs=[('old_contig', ['OLD_00001'])], fna_contigs=['old_contig'])
            result = check_genome_v3(tmp, GENOME_ID)
            self.assertIn('NOT v3', result.summary(GENOME_ID))


class TestCheckV3Deep(TestCase):
    def test_deep_catches_bad_annotation(self):
        with tempfile.TemporaryDirectory() as tmp:
            _setup_genome(tmp,
                          gbk_contigs=_v3_gbk_contigs(),
                          fna_contigs=_v3_fna_contigs(),
                          annotations=[('OLD_00001', 'K00001'), ('OLD_00002', 'K00002')])
            result = check_genome_v3(tmp, GENOME_ID, deep=True)
            self.assertFalse(result.is_v3)
            self.assertTrue(any('locus_tags not v3' in issue for issue in result.issues))

    def test_shallow_ignores_bad_annotation(self):
        with tempfile.TemporaryDirectory() as tmp:
            _setup_genome(tmp,
                          gbk_contigs=_v3_gbk_contigs(),
                          fna_contigs=_v3_fna_contigs(),
                          annotations=[('OLD_00001', 'K00001')])
            result = check_genome_v3(tmp, GENOME_ID, deep=False)
            self.assertTrue(result.is_v3)

    def test_deep_v3_annotation_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            _setup_genome(tmp,
                          gbk_contigs=_v3_gbk_contigs(1, 2),
                          fna_contigs=_v3_fna_contigs(1),
                          annotations=[(f'{GENOME_ID}_000001', 'K00001'),
                                       (f'{GENOME_ID}_000002', 'K00002')])
            result = check_genome_v3(tmp, GENOME_ID, deep=True)
            self.assertTrue(result.is_v3)

    def test_deep_catches_bad_eggnog(self):
        with tempfile.TemporaryDirectory() as tmp:
            _setup_genome(tmp,
                          gbk_contigs=_v3_gbk_contigs(),
                          fna_contigs=_v3_fna_contigs(),
                          eggnog_rows=[('2WKGQ', 'OLD_00001', 'COG0001'),
                                       ('2WKGQ', 'OLD_00002', 'COG0002')])
            result = check_genome_v3(tmp, GENOME_ID, deep=True)
            self.assertFalse(result.is_v3)

    def test_deep_v3_eggnog_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            _setup_genome(tmp,
                          gbk_contigs=_v3_gbk_contigs(1, 2),
                          fna_contigs=_v3_fna_contigs(1),
                          eggnog_rows=[('2WKGQ', f'{GENOME_ID}_000001', 'COG0001'),
                                       ('2WKGQ', f'{GENOME_ID}_000002', 'COG0002')])
            result = check_genome_v3(tmp, GENOME_ID, deep=True)
            self.assertTrue(result.is_v3)


class TestCheckV3PendingFiles(TestCase):
    def test_detects_pending_gbk_v3(self):
        with tempfile.TemporaryDirectory() as tmp:
            _setup_genome(tmp,
                          gbk_contigs=_v3_gbk_contigs(),
                          fna_contigs=_v3_fna_contigs())
            # Simulate a partial upgrade
            open(os.path.join(tmp, f'{GENOME_ID}.gbk.v3'), 'w').close()
            result = check_genome_v3(tmp, GENOME_ID)
            self.assertTrue(result.has_pending_v3_files)
            self.assertTrue(any(f.endswith('.gbk.v3') for f in result.pending_files))
            self.assertIn('PARTIAL UPGRADE', result.summary(GENOME_ID))

    def test_detects_pending_fna_v3(self):
        with tempfile.TemporaryDirectory() as tmp:
            _setup_genome(tmp,
                          gbk_contigs=_v3_gbk_contigs(),
                          fna_contigs=_v3_fna_contigs())
            open(os.path.join(tmp, f'{GENOME_ID}.fna.v3'), 'w').close()
            result = check_genome_v3(tmp, GENOME_ID)
            self.assertTrue(result.has_pending_v3_files)

    def test_deep_detects_pending_annotation_v3(self):
        with tempfile.TemporaryDirectory() as tmp:
            _setup_genome(tmp,
                          gbk_contigs=_v3_gbk_contigs(),
                          fna_contigs=_v3_fna_contigs(),
                          annotations=[(f'{GENOME_ID}_000001', 'K00001')])
            open(os.path.join(tmp, f'{GENOME_ID}.KG.v3'), 'w').close()
            # shallow: not detected (annotations not scanned)
            result_shallow = check_genome_v3(tmp, GENOME_ID, deep=False)
            self.assertFalse(result_shallow.has_pending_v3_files)
            # deep: detected
            result_deep = check_genome_v3(tmp, GENOME_ID, deep=True)
            self.assertTrue(result_deep.has_pending_v3_files)

    def test_no_pending_files_when_clean(self):
        with tempfile.TemporaryDirectory() as tmp:
            _setup_genome(tmp, gbk_contigs=_v3_gbk_contigs(), fna_contigs=_v3_fna_contigs())
            result = check_genome_v3(tmp, GENOME_ID)
            self.assertFalse(result.has_pending_v3_files)
            self.assertEqual(result.pending_files, [])


class TestCheckV3SubdirectoryLayout(TestCase):
    """check_genome_v3 should work regardless of whether files are in subdirs."""

    def test_subdirectory_layout(self):
        with tempfile.TemporaryDirectory() as tmp:
            cds_dir = os.path.join(tmp, '2_cds')
            asm_dir = os.path.join(tmp, '1_assembly')
            os.makedirs(cds_dir)
            os.makedirs(asm_dir)

            gbk_path = os.path.join(cds_dir, f'{GENOME_ID}.gbk')
            fna_path = os.path.join(asm_dir, f'{GENOME_ID}.fna')

            _write_gbk(gbk_path, _v3_gbk_contigs())
            _write_fna(fna_path, _v3_fna_contigs())

            genome_json = {
                'identifier': GENOME_ID,
                'cds_tool_gbk_file': f'2_cds/{GENOME_ID}.gbk',
                'assembly_fasta_file': f'1_assembly/{GENOME_ID}.fna',
            }
            with open(os.path.join(tmp, 'genome.json'), 'w') as f:
                json.dump(genome_json, f)

            result = check_genome_v3(tmp, GENOME_ID)
            self.assertTrue(result.is_v3)
