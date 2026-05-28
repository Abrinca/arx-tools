import os
import tempfile
import unittest

from arx_tools.rename_fasta import FastaFile

ROOT = os.path.dirname(os.path.dirname(__file__))

SAMPLE_FAA = """\
>OLD_000001 hypothetical protein
MAST
>OLD_000002 transposase
MKVL
"""

LT_MAP = {
    'OLD_000001': 'NEW_000001',
    'OLD_000002': 'NEW_000002',
}

fastas = [
    f'{ROOT}/test-data/prokka-bad/PROKKA_08112021.',
    f'{ROOT}/test-data/prokka-good/PROKKA_08112021.',
    f'{ROOT}/test-data/pgap-bad/annot.',
    f'{ROOT}/test-data/pgap-good/annot.'
]
fastas = [f + suffix for suffix in ['faa', 'ffn'] for f in fastas][:-2]


class Test(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def test_parse_fasta_header(self):
        for fasta in fastas:
            with open(fasta) as f:
                locus_tag_prefix, gene_id = FastaFile.parse_fasta_header(f.readline())
            self.assertIn(member=locus_tag_prefix, container=['tmp_', 'STRAIN.1_'])
            self.assertIn(member=gene_id, container=['00001', '000289'])

    def test_detect_locus_tag_prefix(self):
        for fasta in fastas:
            locus_tag_prefix = FastaFile(fasta).detect_locus_tag_prefix()
            self.assertIn(member=locus_tag_prefix, container=['tmp_', 'STRAIN.1_'])

    def test_rename(self):
        for fasta in fastas:
            out = os.path.join(self.tmp, '_'.join(fasta.split(os.sep)[-2:]))
            FastaFile(fasta).rename(new_locus_tag_prefix='YOLO_', out=out, validate=True)
            with open(fasta) as f_old, open(out) as f_new:
                content_old = f_old.read()
                content_new = f_new.read()
            self.assertNotIn(member='tmp', container=content_new)
            self.assertNotIn(member='STRAIN.1', container=content_new)
            self.assertEqual(
                first=content_new.count('YOLO_'),
                second=max(content_old.count('tmp_'), content_old.count('STRAIN.1_'))
            )

    def test_rename_by_map(self):
        src = os.path.join(self.tmp, 'input.faa')
        out = os.path.join(self.tmp, 'renamed.faa')
        with open(src, 'w') as f:
            f.write(SAMPLE_FAA)
        FastaFile(src).rename_by_map(out=out, lt_map=LT_MAP, update_path=False)
        with open(out) as f:
            content = f.read()
        self.assertNotIn('OLD_000001', content)
        self.assertNotIn('OLD_000002', content)
        self.assertIn('>NEW_000001 hypothetical protein\n', content)
        self.assertIn('>NEW_000002 transposase\n', content)
        self.assertIn('MAST\n', content)
        self.assertIn('MKVL\n', content)
