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

# Prokka GFF with gnl|X| contig IDs (as produced when input FASTA uses NCBI-style headers)
SAMPLE_GFF_GNL = """\
##gff-version 3
##sequence-region gnl|C|BARE_1 1 40066
##sequence-region gnl|C|BARE_2 1 39909
gnl|C|BARE_1\tprokka\tgene\t32\t637\t.\t-\t.\tID=OLD_000001_gene;locus_tag=OLD_000001
gnl|C|BARE_1\tProdigal:002006\tCDS\t32\t637\t.\t-\t0\tID=OLD_000001;Parent=OLD_000001_gene;locus_tag=OLD_000001;product=hypothetical protein
gnl|C|BARE_2\tprokka\tgene\t100\t500\t.\t+\t.\tID=OLD_000002_gene;locus_tag=OLD_000002
gnl|C|BARE_2\tProdigal:002006\tCDS\t100\t500\t.\t+\t0\tID=OLD_000002;Parent=OLD_000002_gene;locus_tag=OLD_000002;product=hypothetical protein
"""

SAMPLE_GFF_GNL_WITH_FASTA = SAMPLE_GFF_GNL + """\
##FASTA
>gnl|C|BARE_1
ATGC
>gnl|C|BARE_2
TTTT
"""

LT_MAP = {
    'OLD_000001': 'NEW_000001',
    'OLD_000002': 'NEW_000002',
}

CONTIG_ID_MAP = {
    'BARE_1': 'GENOME_scf1',
    'BARE_2': 'GENOME_scf2',
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

    def test_rename_by_map_with_gnl_contig_ids(self):
        """GFFs from Prokka on NCBI-style input use gnl|X|bare_id contig IDs.
        rename_by_map must update ##sequence-region headers and column 0."""
        src = os.path.join(self.tmp, 'input_gnl.gff')
        out = os.path.join(self.tmp, 'renamed_gnl.gff')
        with open(src, 'w') as f:
            f.write(SAMPLE_GFF_GNL)
        GffFile(src).rename_by_map(out=out, lt_map=LT_MAP, contig_id_map=CONTIG_ID_MAP, update_path=False)
        with open(out) as f:
            content = f.read()
        # locus tags renamed
        self.assertIn('locus_tag=NEW_000001', content)
        self.assertIn('locus_tag=NEW_000002', content)
        # contig IDs in data column 0 renamed
        self.assertIn('GENOME_scf1\t', content)
        self.assertIn('GENOME_scf2\t', content)
        self.assertNotIn('gnl|C|BARE_1', content)
        self.assertNotIn('gnl|C|BARE_2', content)
        # ##sequence-region headers renamed
        self.assertIn('##sequence-region GENOME_scf1 ', content)
        self.assertIn('##sequence-region GENOME_scf2 ', content)
        # original lengths preserved
        self.assertIn('##sequence-region GENOME_scf1 1 40066', content)
        self.assertIn('##sequence-region GENOME_scf2 1 39909', content)

    def test_rename_by_map_gnl_embedded_fasta_renamed(self):
        """Embedded ##FASTA contig headers must also be renamed when contig_id_map is set."""
        src = os.path.join(self.tmp, 'input_gnl_fasta.gff')
        out = os.path.join(self.tmp, 'renamed_gnl_fasta.gff')
        with open(src, 'w') as f:
            f.write(SAMPLE_GFF_GNL_WITH_FASTA)
        GffFile(src).rename_by_map(out=out, lt_map=LT_MAP, contig_id_map=CONTIG_ID_MAP, update_path=False)
        with open(out) as f:
            content = f.read()
        self.assertNotIn('gnl|C|BARE_1', content)
        self.assertNotIn('gnl|C|BARE_2', content)
        self.assertIn('>GENOME_scf1\n', content)
        self.assertIn('>GENOME_scf2\n', content)
        self.assertIn('ATGC\n', content)
        self.assertIn('TTTT\n', content)
