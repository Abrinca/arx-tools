import tempfile
from unittest import TestCase

import os
from arx_tools.rename_custom_annotations import *

CA_CONTENT = """\
OLD_00001\tK00001\t1.0
OLD_00002\tK00002\t0.9
OLD_00003\tK00003\t0.8
"""

ROOT = os.path.dirname(os.path.dirname(__file__))
TMPFILE = '/tmp/renamed_custom_annotations.KG'

custom_files = [
    f'{ROOT}/test-data/prokka-bad/custom_annotations.KG',
]


def cleanup():
    if os.path.isfile(TMPFILE):
        os.remove(TMPFILE)


class Test(TestCase):
    def test_detect_locus_tag_prefix(self):
        for custom_file in custom_files:
            locus_tag_prefix = CustomAnnotationFile(custom_file).detect_locus_tag_prefix()
            self.assertEqual(locus_tag_prefix, 'tmp_')

    def test_validate_locus_tags(self):
        for custom_file in custom_files:
            CustomAnnotationFile(custom_file).validate_locus_tags(locus_tag_prefix=None)
            CustomAnnotationFile(custom_file).validate_locus_tags(locus_tag_prefix='tmp_')
            with self.assertRaises(AssertionError):
                CustomAnnotationFile(custom_file).validate_locus_tags(locus_tag_prefix='xxx_')

    def test_rename_custom_annotations(self):
        for custom_file in custom_files:
            cleanup()
            CustomAnnotationFile(custom_file).rename(new_locus_tag_prefix='YOLO_', out=TMPFILE, validate=True)
            with open(TMPFILE) as f:
                content = f.read()
            count = content.count('YOLO_')
            self.assertNotIn(member='tmp', container=content)
            self.assertEqual(count, 3)

    def test_rename_by_map(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.KG', delete=False) as f:
            f.write(CA_CONTENT)
            infile = f.name
        outfile = infile + '.out'
        try:
            lt_map = {'OLD_00001': 'NEW_00001', 'OLD_00002': 'NEW_00002', 'OLD_00003': 'NEW_00003'}
            CustomAnnotationFile(infile).rename_by_map(out=outfile, lt_map=lt_map, update_path=False)
            with open(outfile) as f:
                result = f.read()
            self.assertIn('NEW_00001', result)
            self.assertIn('NEW_00002', result)
            self.assertIn('NEW_00003', result)
            self.assertNotIn('OLD_', result)
            self.assertIn('K00001', result)  # rest of line preserved
        finally:
            os.remove(infile)
            if os.path.isfile(outfile):
                os.remove(outfile)

    def test_rename_by_map_missing_key(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.KG', delete=False) as f:
            f.write(CA_CONTENT)
            infile = f.name
        outfile = infile + '.out'
        try:
            lt_map = {'OLD_00001': 'NEW_00001'}  # missing OLD_00002 and OLD_00003
            with self.assertRaises(AssertionError):
                CustomAnnotationFile(infile).rename_by_map(out=outfile, lt_map=lt_map)
        finally:
            os.remove(infile)
            if os.path.isfile(outfile):
                os.remove(outfile)

    @classmethod
    def tearDownClass(cls) -> None:
        cleanup()
