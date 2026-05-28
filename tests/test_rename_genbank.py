import tempfile
from unittest import TestCase

import os
from arx_tools.rename_genbank import *
from arx_tools.rename_fasta import FastaFile

# Re-import as classes after the wildcard (which pulls in Bio.SeqRecord etc. as modules)
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO


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


def _write_fna(path: str, contigs: list[tuple[str, str]]) -> None:
    """Write a minimal FASTA file. contigs: [(contig_id, description), ...]"""
    with open(path, 'w') as f:
        for contig_id, desc in contigs:
            header = f'>{contig_id} {desc}' if desc else f'>{contig_id}'
            f.write(f'{header}\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n')

ROOT = os.path.dirname(os.path.dirname(__file__))
TMPFILE = '/tmp/renamed_gbk.gbk'

gbks = [
    f'{ROOT}/test-data/prokka-bad/PROKKA_08112021.gbk',
    f'{ROOT}/test-data/prokka-good/PROKKA_08112021.gbk',
    f'{ROOT}/test-data/pgap-bad/annot.gbk',
    f'{ROOT}/test-data/pgap-good/annot.gbk'
]


def cleanup():
    if os.path.isfile(TMPFILE):
        os.remove(TMPFILE)


class Test(TestCase):
    def test_detect_locus_tag_prefix(self):
        for gbk in gbks:
            strain, locus_tag_prefix = GenBankFile(gbk).detect_strain_locus_tag_prefix()
            self.assertIn(member=strain, container=['replaceme', 'STRAIN'])
            self.assertIn(member=locus_tag_prefix, container=['tmp_', 'STRAIN.1_'])

    def test_get_taxid(self):
        for gbk in gbks:
            self.assertEqual(GenBankFile(gbk).taxid(), 2097)

    def test_rename(self):
        for gbk in gbks:
            cleanup()
            GenBankFile(gbk).rename(new_locus_tag_prefix='YOLO_', out=TMPFILE, validate=True)
            with open(TMPFILE) as f:
                content = f.read()
            count = content.count('YOLO_')
            self.assertNotIn(member='tmp_', container=content)
            self.assertNotIn(member='STRAIN.1_', container=content)
            self.assertGreater(a=count, b=1000)

    def test_rename_reindex(self):
        for gbk in gbks:
            cleanup()
            GenBankFile(gbk).rename(
                new_locus_tag_prefix='YOLO_', out=TMPFILE, validate=True,
                scf_prefix='XX_scf', scf_leading_zeroes=10
            )
            with open(TMPFILE) as f:
                content = f.read()
            count = content.count('XX_scf0')
            self.assertGreater(a=count, b=0)

    def test_get_metadata(self):
        for gbk in gbks:
            organism_data, genome_data = GenBankFile(gbk).metadata()
            self.assertIn(genome_data['cds_tool'], container=['PGAP', 'prokka'])
            self.assertIn(genome_data['cds_tool_date'], container=['2021-08-10', '2021-08-11'])
            self.assertIn(genome_data['cds_tool_version'], container=['2021-07-01.build5508', '1.14.5'])

    def test_get_contig_ids(self):
        with tempfile.NamedTemporaryFile(suffix='.gbk', delete=False) as f:
            path = f.name
        try:
            _write_gbk(path, [('contig_B', ['B_00001']), ('contig_A', ['A_00001'])])
            ids = GenBankFile(path).get_contig_ids()
            self.assertEqual(ids, ['contig_B', 'contig_A'])
        finally:
            os.remove(path)

    def test_normalize_updates_protein_id(self):
        """protein_id=C:{old_lt} must be updated to C:{new_lt} alongside locus_tag."""
        with (tempfile.NamedTemporaryFile(suffix='.gbk', delete=False) as in_f,
              tempfile.NamedTemporaryFile(suffix='.gbk', delete=False) as out_f):
            in_path, out_path = in_f.name, out_f.name
        try:
            rec = SeqRecord(Seq('ATCGATCG' * 10), id='contig1', name='contig1', description='')
            rec.annotations['molecule_type'] = 'DNA'
            cds = SeqFeature(FeatureLocation(0, 9, strand=1), type='CDS')
            cds.qualifiers['locus_tag'] = ['OLD_00001']
            cds.qualifiers['protein_id'] = ['C:OLD_00001']
            rec.features.append(cds)
            # A second feature without protein_id — must not crash
            gene = SeqFeature(FeatureLocation(0, 9, strand=1), type='gene')
            gene.qualifiers['locus_tag'] = ['OLD_00001']
            rec.features.append(gene)
            with open(in_path, 'w') as f:
                SeqIO.write([rec], f, 'genbank')

            GenBankFile(in_path).normalize(out=out_path, genome_id='GENOME')

            with open(out_path) as f:
                result = list(SeqIO.parse(f, 'genbank'))
            feat = next(ft for ft in result[0].features if ft.type == 'CDS')
            self.assertEqual(feat.qualifiers['locus_tag'], ['GENOME_000001'])
            self.assertEqual(feat.qualifiers['protein_id'], ['C:GENOME_000001'])
        finally:
            for p in (in_path, out_path):
                if os.path.isfile(p):
                    os.remove(p)

    def test_contig_order_fna_gbk_same_contigs_different_order(self):
        """Canonical IDs follow GBK order; FNA contigs are renamed to match regardless of FNA order."""
        with (tempfile.NamedTemporaryFile(suffix='.gbk', delete=False) as gbk_f,
              tempfile.NamedTemporaryFile(suffix='.fna', delete=False) as fna_f,
              tempfile.NamedTemporaryFile(suffix='.gbk', delete=False) as out_gbk_f,
              tempfile.NamedTemporaryFile(suffix='.fna', delete=False) as out_fna_f):
            gbk_path, fna_path = gbk_f.name, fna_f.name
            out_gbk, out_fna = out_gbk_f.name, out_fna_f.name

        try:
            # GBK order: [contig_B, contig_A]  →  canonical: B→scf00001, A→scf00002
            _write_gbk(gbk_path, [
                ('contig_B', ['OLD_00001']),
                ('contig_A', ['OLD_00002']),
            ])
            # FNA order: [contig_A, contig_B] deliberately reversed
            _write_fna(fna_path, [
                ('contig_A', 'topology=linear coverage=50x'),
                ('contig_B', 'topology=circular coverage=200x'),
            ])

            gbk = GenBankFile(gbk_path)
            fna = FastaFile(fna_path)
            genome = 'GENOME'
            contig_format = '_scf{n:05d}'

            gbk_contig_ids = gbk.get_contig_ids()
            fna_contig_ids = fna.get_contig_ids()
            self.assertEqual(set(gbk_contig_ids), set(fna_contig_ids))

            contig_id_map = {
                gbk_id: f'{genome}{contig_format.format(n=i + 1)}'
                for i, gbk_id in enumerate(gbk_contig_ids)
            }

            # rename FNA in FNA order
            fna.rename_contig_ids(out=out_fna, new_ids=[contig_id_map[i] for i in fna_contig_ids])
            # normalize GBK in GBK order
            canonical_ids = [contig_id_map[i] for i in gbk_contig_ids]
            gbk.normalize(out=out_gbk, genome_id=genome, contig_ids=canonical_ids, contig_format=contig_format)

            renamed_fna_ids = FastaFile(out_fna).get_contig_ids()
            renamed_gbk_ids = GenBankFile(out_gbk).get_contig_ids()

            # contig_A is 2nd in GBK → scf00002; contig_B is 1st in GBK → scf00001
            self.assertEqual(renamed_gbk_ids, ['GENOME_scf00001', 'GENOME_scf00002'])
            # FNA was [contig_A, contig_B] → [scf00002, scf00001]
            self.assertEqual(renamed_fna_ids, ['GENOME_scf00002', 'GENOME_scf00001'])

            # metadata preserved in FNA headers
            with open(out_fna) as f:
                fna_content = f.read()
            self.assertIn('topology=linear coverage=50x', fna_content)
            self.assertIn('topology=circular coverage=200x', fna_content)
        finally:
            for p in (gbk_path, fna_path, out_gbk, out_fna):
                if os.path.isfile(p):
                    os.remove(p)

    def test_contig_order_mismatch_raises(self):
        """AssertionError when FNA and GBK contain different contig IDs."""
        with (tempfile.NamedTemporaryFile(suffix='.gbk', delete=False) as gbk_f,
              tempfile.NamedTemporaryFile(suffix='.fna', delete=False) as fna_f):
            gbk_path, fna_path = gbk_f.name, fna_f.name
        try:
            _write_gbk(gbk_path, [('contig_A', ['OLD_00001']), ('contig_B', ['OLD_00002'])])
            _write_fna(fna_path, [('contig_A', ''), ('contig_X', '')])  # contig_X not in GBK

            gbk_ids = set(GenBankFile(gbk_path).get_contig_ids())
            fna_ids = set(FastaFile(fna_path).get_contig_ids())
            self.assertNotEqual(gbk_ids, fna_ids)
            with self.assertRaises(AssertionError):
                assert gbk_ids == fna_ids, (
                    f'FNA and GBK contain different contig IDs.\n'
                    f'  FNA only: {fna_ids - gbk_ids}\n'
                    f'  GBK only: {gbk_ids - fna_ids}'
                )
        finally:
            for p in (gbk_path, fna_path):
                if os.path.isfile(p):
                    os.remove(p)

    @classmethod
    def tearDownClass(cls) -> None:
        cleanup()
