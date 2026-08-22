from unittest import TestCase

import os
import logging
import tempfile
from arx_tools.genbank_to_fasta import GenBankToFasta

logging.basicConfig(level=logging.INFO)

ROOT = os.path.dirname(os.path.dirname(__file__))

GENBANK_FILES = [
    (f'{ROOT}/test-data/pgap-bad/annot.gbk', 'tmp_'),
    (f'{ROOT}/test-data/prokka-bad/PROKKA_08112021.gbk', 'tmp_'),
    (f'{ROOT}/test-data/ncbi-download/GCF_005864195.1.gbk', 'FEZ40_RS')
]

FFN_CHARS = set('NATGC')
FAA_CHARS = set('YPMLWDKHEARGTQSNVFCI')

TMPFILE = '/tmp/genbank_to_fasta.fasta'


def cleanup():
    if os.path.isfile(TMPFILE):
        os.remove(TMPFILE)


class TestGenBankToFasta(TestCase):
    def setUp(self) -> None:
        for gbk, locus_tag_prefix in GENBANK_FILES:
            assert os.path.isfile(gbk), f'File does not exist: {gbk}'

    def convert_tester(self, gbk, format):
        cleanup()
        GenBankToFasta.convert(gbk=gbk, out=TMPFILE, format=format, strict=True)
        cleanup()

    def test_convert(self):
        for gbk, locus_tag_prefix in GENBANK_FILES:
            for format in ['faa', 'ffn']:
                self.convert_tester(gbk=gbk, format=format)

    def long_fasta_tester(self, gbk, format, locus_tag_prefix):
        ALLOWED_CHARS = FAA_CHARS if format == 'faa' else FFN_CHARS
        for entry in GenBankToFasta._long_fasta_generator(gbk=gbk, format=format, strict=True):
            self.assertEqual(entry.count('\n'), 1, f'Expect 1 newlines! {entry=}')
            header, sequence = entry.split('\n', 1)
            self.assertTrue(header.startswith(f'>{locus_tag_prefix}'), f'Header does not start with >{locus_tag_prefix}: {entry=}')
            self.assertTrue(set(sequence).issubset(ALLOWED_CHARS),
                            msg=f'FASTA sequences may only contain {ALLOWED_CHARS}. {set(sequence)=} {sequence=}')

    def test_long_fasta_generator(self):
        for gbk, locus_tag_prefix in GENBANK_FILES:
            for format in ['faa', 'ffn']:
                logging.info(f'testing {format=} {gbk=}')
                self.long_fasta_tester(gbk=gbk, format=format, locus_tag_prefix=locus_tag_prefix)

    def short_fasta_tester(self, gbk, format, locus_tag_prefix):
        ALLOWED_CHARS = FAA_CHARS if format == 'faa' else FFN_CHARS
        short_fasta_generator = GenBankToFasta._short_fasta_generator(gbk=gbk, format=format, strict=True)

        for entry in short_fasta_generator:
            self.assertGreater(entry.count('\n'), 0, f'FASTA entry must at least contain 1 newline! {entry=}')
            header, *sequence = entry.split('\n')

            self.assertTrue(header.startswith(f'>{locus_tag_prefix}'), f'Header does not start with >{locus_tag_prefix}: {entry=}')

            for line in sequence[:-3]:
                self.assertEqual(len(line), 60, msg=f'Short FASTA lines must be 60 chars long. {line=}')
                self.assertTrue(set(line).issubset(ALLOWED_CHARS), msg=f'FASTA sequences may only contain {ALLOWED_CHARS}. {set(line)=} {entry=}')

            last_line = sequence[-2]
            self.assertTrue(set(last_line).issubset(ALLOWED_CHARS),
                            msg=f'FASTA sequences may only contain {ALLOWED_CHARS}. {set(last_line)=} {entry=}')
            self.assertTrue(0 < len(last_line) <= 60, msg=f'last FASTA line must be between 1 and 60 chars long. {last_line=}')

            self.assertTrue(sequence[-1] == '', f'FASTA entries must end with a newline. {entry=}')

    def test_short_fasta_generator(self):
        for gbk, locus_tag_prefix in GENBANK_FILES:
            for format in ['faa', 'ffn']:
                logging.info(f'testing {format=} {gbk=}')
                self.short_fasta_tester(gbk=gbk, format=format, locus_tag_prefix=locus_tag_prefix)



class TestCreateGff(TestCase):
    def test_create_gff(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation

        locus_tags = ['MYORG_000001', 'MYORG_000002']
        rec = SeqRecord(Seq('ATCGATCGATCG' * 10), id='MYORG_scf1', name='MYORG_scf1', description='')
        rec.annotations['molecule_type'] = 'DNA'
        for i, lt in enumerate(locus_tags):
            gene = SeqFeature(FeatureLocation(i * 30, i * 30 + 30, strand=1), type='gene')
            gene.qualifiers['locus_tag'] = [lt]
            cds = SeqFeature(FeatureLocation(i * 30, i * 30 + 30, strand=1), type='CDS')
            cds.qualifiers['locus_tag'] = [lt]
            cds.qualifiers['product'] = [f'hypothetical protein {lt}']
            rec.features += [gene, cds]

        with (tempfile.NamedTemporaryFile(suffix='.gbk', mode='w', delete=False) as gbk_f,
              tempfile.NamedTemporaryFile(suffix='.gff', delete=False) as gff_f):
            gbk_path, gff_path = gbk_f.name, gff_f.name

        try:
            from Bio import SeqIO
            with open(gbk_path, 'w') as f:
                SeqIO.write([rec], f, 'genbank')

            GenBankToFasta.create_gff(gbk=gbk_path, out=gff_path)

            with open(gff_path) as f:
                content = f.read()

            self.assertIn('##gff-version 3', content)
            self.assertIn('MYORG_000001', content)
            self.assertIn('MYORG_000002', content)

            feature_lines = [l for l in content.splitlines()
                             if l and not l.startswith('#')]
            self.assertGreater(len(feature_lines), 0)
            for line in feature_lines:
                cols = line.split('\t')
                self.assertEqual(len(cols), 9, f'Expected 9 GFF columns: {line!r}')

            # CDS features must have a Parent pointing to the gene
            cds_lines = [l for l in feature_lines if l.split('\t')[2] == 'CDS']
            self.assertGreater(len(cds_lines), 0)
            for line in cds_lines:
                self.assertIn('Parent=', line, f'CDS line missing Parent: {line!r}')
                self.assertIn('locus_tag=', line)

            # gene features must have ID attributes
            gene_lines = [l for l in feature_lines if l.split('\t')[2] == 'gene']
            self.assertGreater(len(gene_lines), 0)
            for line in gene_lines:
                self.assertIn('ID=', line, f'gene line missing ID: {line!r}')
        finally:
            for p in (gbk_path, gff_path):
                if os.path.exists(p):
                    os.unlink(p)
