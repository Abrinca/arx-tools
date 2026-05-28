import os
import tempfile
import unittest

from arx_tools.rename_gff import GffFile

ROOT = os.path.dirname(os.path.dirname(__file__))

gffs = [
    f'{ROOT}/test-data/prokka-bad/PROKKA_08112021.gff',
    f'{ROOT}/test-data/prokka-good/PROKKA_08112021.gff',
    f'{ROOT}/test-data/pgap-bad/annot.gff',
    f'{ROOT}/test-data/pgap-good/annot.gff'
]

SAMPLE_GFF = """\
##gff-version 3
##sequence-region GENOME_scf1 1 5000
GENOME_scf1\tProkka\tgene\t1\t900\t.\t+\t.\tID=gene-OLD_000001;Name=OLD_000001;locus_tag=OLD_000001
GENOME_scf1\tProkka\tCDS\t1\t900\t.\t+\t0\tID=cds-OLD_000001;Parent=gene-OLD_000001;locus_tag=OLD_000001;product=hypothetical protein
GENOME_scf1\tProkka\tgene\t1000\t1900\t.\t+\t.\tID=gene-OLD_000002;Name=OLD_000002;locus_tag=OLD_000002
GENOME_scf1\tProkka\tCDS\t1000\t1900\t.\t+\t0\tID=cds-OLD_000002;Parent=gene-OLD_000002;locus_tag=OLD_000002;product=hypothetical protein
##FASTA
>GENOME_scf1
ATGC
"""

LT_MAP = {
    'OLD_000001': 'NEW_000001',
    'OLD_000002': 'NEW_000002',
}


class Test(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def test_detect_locus_tag_prefix(self):
        for gff in gffs:
            locus_tag_prefix = GffFile(gff).detect_locus_tag_prefix()
            self.assertIn(member=locus_tag_prefix, container=['tmp_', 'STRAIN.1_'])

    def test_rename(self):
        for gff in gffs:
            out = os.path.join(self.tmp, '_'.join(gff.split(os.sep)[-2:]))
            GffFile(gff).rename(new_locus_tag_prefix='YOLO_', out=out, validate=True)
            with open(out) as f:
                content = f.read()
            count = content.count('YOLO_')
            self.assertNotIn(member='tmp', container=content, msg=f'Found "tmp" in renamed {gff=}!')
            self.assertNotIn(member='STRAIN.1', container=content, msg=f'Found "STRAIN.1" in renamed {gff=}!')
            self.assertGreater(a=count, b=1000)

    def test_rename_by_map(self):
        src = os.path.join(self.tmp, 'input.gff')
        out = os.path.join(self.tmp, 'renamed.gff')
        with open(src, 'w') as f:
            f.write(SAMPLE_GFF)
        GffFile(src).rename_by_map(out=out, lt_map=LT_MAP, update_path=False)
        with open(out) as f:
            content = f.read()
        self.assertNotIn('OLD_000001', content)
        self.assertNotIn('OLD_000002', content)
        self.assertIn('locus_tag=NEW_000001', content)
        self.assertIn('locus_tag=NEW_000002', content)
        self.assertIn('ID=gene-NEW_000001', content)
        self.assertIn('Parent=gene-NEW_000001', content)
        self.assertIn('##FASTA\n', content)
        self.assertIn('>GENOME_scf1\n', content)
        self.assertIn('##gff-version 3\n', content)
