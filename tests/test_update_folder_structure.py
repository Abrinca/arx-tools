import json
import os
import tarfile
import tempfile
from unittest import TestCase

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

from arx_tools.update_folder_structure import (
    _apply_lt_map_to_file,
    _apply_contig_map_to_fna,
    _promote_v3_files,
)
from arx_tools.check_v3 import check_genome_v3
from arx_tools.rename_genbank import GenBankFile

GENOME_ID = 'FAM1079-i1-1.1'


def _write_gbk(path: str, contigs: list[tuple[str, list[str]]]) -> None:
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
    with open(path, 'w') as f:
        for contig_id in contigs:
            f.write(f'>{contig_id} some description\nATCGATCGATCG\n')


def _write_annotation(path: str, rows: list[tuple]) -> None:
    with open(path, 'w') as f:
        f.write('# comment\n')
        for row in rows:
            f.write('\t'.join(row) + '\n')


class TestApplyLtMapToFile(TestCase):
    def test_basic_rename(self):
        lt_map = {'OLD_000001': 'NEW_000001', 'OLD_000002': 'NEW_000002'}
        with tempfile.NamedTemporaryFile(mode='w', suffix='.KG', delete=False) as src:
            src.write('# comment\nOLD_000001\tK00001\nOLD_000002\tK00002\n')
            src_path = src.name
        dst_path = src_path + '.out'
        try:
            _apply_lt_map_to_file(src_path, dst_path, lt_map)
            with open(dst_path) as f:
                content = f.read()
            self.assertIn('NEW_000001\tK00001', content)
            self.assertIn('NEW_000002\tK00002', content)
            self.assertNotIn('OLD_', content)
            self.assertIn('# comment', content)
        finally:
            os.unlink(src_path)
            if os.path.exists(dst_path):
                os.unlink(dst_path)

    def test_eggnog_prefix_preserved(self):
        lt_map = {'OLD_000001': 'NEW_000001'}
        with tempfile.NamedTemporaryFile(mode='w', suffix='.annotations', delete=False) as src:
            src.write('# eggnog\n2WKGQ|OLD_000001\tCOG0001\n')
            src_path = src.name
        dst_path = src_path + '.out'
        try:
            _apply_lt_map_to_file(src_path, dst_path, lt_map, is_eggnog=True)
            with open(dst_path) as f:
                content = f.read()
            self.assertIn('2WKGQ|NEW_000001', content)
            self.assertNotIn('OLD_', content)
        finally:
            os.unlink(src_path)
            if os.path.exists(dst_path):
                os.unlink(dst_path)

    def test_unmapped_tag_preserved(self):
        """Tags not in lt_map are left unchanged."""
        lt_map = {'OLD_000001': 'NEW_000001'}
        with tempfile.NamedTemporaryFile(mode='w', suffix='.KG', delete=False) as src:
            src.write('OLD_000001\tK00001\nUNKNOWN_TAG\tK00002\n')
            src_path = src.name
        dst_path = src_path + '.out'
        try:
            _apply_lt_map_to_file(src_path, dst_path, lt_map)
            with open(dst_path) as f:
                content = f.read()
            self.assertIn('NEW_000001', content)
            self.assertIn('UNKNOWN_TAG', content)
        finally:
            os.unlink(src_path)
            if os.path.exists(dst_path):
                os.unlink(dst_path)

    def test_source_file_unchanged(self):
        """src is never modified — only dst is written."""
        lt_map = {'OLD_000001': 'NEW_000001'}
        original = 'OLD_000001\tK00001\n'
        with tempfile.NamedTemporaryFile(mode='w', suffix='.KG', delete=False) as src:
            src.write(original)
            src_path = src.name
        dst_path = src_path + '.out'
        try:
            _apply_lt_map_to_file(src_path, dst_path, lt_map)
            with open(src_path) as f:
                self.assertEqual(f.read(), original)
        finally:
            os.unlink(src_path)
            if os.path.exists(dst_path):
                os.unlink(dst_path)


class TestApplyContigMapToFna(TestCase):
    def test_renames_headers(self):
        contig_map = {'old_scf1': 'NEW_scf1', 'old_scf2': 'NEW_scf2'}
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fna', delete=False) as src:
            src.write('>old_scf1 desc1\nATCG\n>old_scf2 desc2\nATCG\n')
            src_path = src.name
        dst_path = src_path + '.out'
        try:
            count = _apply_contig_map_to_fna(src_path, dst_path, contig_map)
            self.assertEqual(count, 2)
            with open(dst_path) as f:
                content = f.read()
            self.assertIn('>NEW_scf1 desc1', content)
            self.assertIn('>NEW_scf2 desc2', content)
            self.assertNotIn('>old_scf', content)
        finally:
            os.unlink(src_path)
            if os.path.exists(dst_path):
                os.unlink(dst_path)

    def test_description_preserved(self):
        contig_map = {'old_scf1': 'NEW_scf1'}
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fna', delete=False) as src:
            src.write('>old_scf1 topology=linear coverage=50x\nATCG\n')
            src_path = src.name
        dst_path = src_path + '.out'
        try:
            _apply_contig_map_to_fna(src_path, dst_path, contig_map)
            with open(dst_path) as f:
                content = f.read()
            self.assertIn('topology=linear coverage=50x', content)
        finally:
            os.unlink(src_path)
            if os.path.exists(dst_path):
                os.unlink(dst_path)

    def test_unmapped_contig_unchanged(self):
        contig_map = {'old_scf1': 'NEW_scf1'}
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fna', delete=False) as src:
            src.write('>old_scf1\nATCG\n>unlisted_scf\nATCG\n')
            src_path = src.name
        dst_path = src_path + '.out'
        try:
            count = _apply_contig_map_to_fna(src_path, dst_path, contig_map)
            self.assertEqual(count, 1)
            with open(dst_path) as f:
                content = f.read()
            self.assertIn('>unlisted_scf', content)
        finally:
            os.unlink(src_path)
            if os.path.exists(dst_path):
                os.unlink(dst_path)

    def test_source_file_unchanged(self):
        contig_map = {'old_scf1': 'NEW_scf1'}
        original = '>old_scf1\nATCG\n'
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fna', delete=False) as src:
            src.write(original)
            src_path = src.name
        dst_path = src_path + '.out'
        try:
            _apply_contig_map_to_fna(src_path, dst_path, contig_map)
            with open(src_path) as f:
                self.assertEqual(f.read(), original)
        finally:
            os.unlink(src_path)
            if os.path.exists(dst_path):
                os.unlink(dst_path)


class TestPromoteV3Files(TestCase):
    def test_promotes_and_archives(self):
        with tempfile.TemporaryDirectory() as tmp:
            orig = os.path.join(tmp, 'file.gbk')
            v3 = orig + '.v3'
            with open(orig, 'w') as f:
                f.write('original')
            with open(v3, 'w') as f:
                f.write('v3 content')

            archive = _promote_v3_files({v3: orig}, genome_dir=tmp, genome_id='TEST')

            self.assertFalse(os.path.exists(v3))
            self.assertTrue(os.path.exists(orig))
            with open(orig) as f:
                self.assertEqual(f.read(), 'v3 content')
            self.assertTrue(os.path.exists(archive))
            with tarfile.open(archive, 'r:gz') as tar:
                names = tar.getnames()
            self.assertIn('file.gbk', names)

    def test_archive_contains_original_content(self):
        with tempfile.TemporaryDirectory() as tmp:
            orig = os.path.join(tmp, 'file.gbk')
            v3 = orig + '.v3'
            with open(orig, 'w') as f:
                f.write('original content')
            with open(v3, 'w') as f:
                f.write('v3 content')

            archive = _promote_v3_files({v3: orig}, genome_dir=tmp, genome_id='TEST')

            with tarfile.open(archive, 'r:gz') as tar:
                member = tar.extractfile('file.gbk')
                self.assertEqual(member.read().decode(), 'original content')

    def test_promotes_multiple_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            files = {}
            for name in ('a.gbk', 'b.fna', 'c.faa'):
                orig = os.path.join(tmp, name)
                v3 = orig + '.v3'
                with open(orig, 'w') as f:
                    f.write(f'orig_{name}')
                with open(v3, 'w') as f:
                    f.write(f'v3_{name}')
                files[v3] = orig

            archive = _promote_v3_files(files, genome_dir=tmp, genome_id='TEST')

            for v3, orig in files.items():
                self.assertFalse(os.path.exists(v3))
                self.assertTrue(os.path.exists(orig))
            with tarfile.open(archive, 'r:gz') as tar:
                archived = set(tar.getnames())
            self.assertEqual(archived, {'a.gbk', 'b.fna', 'c.faa'})

    def test_subdirectory_paths_in_archive(self):
        """Files in subdirs are stored with relative paths inside the archive."""
        with tempfile.TemporaryDirectory() as tmp:
            subdir = os.path.join(tmp, '2_cds')
            os.makedirs(subdir)
            orig = os.path.join(subdir, 'genome.gbk')
            v3 = orig + '.v3'
            with open(orig, 'w') as f:
                f.write('original')
            with open(v3, 'w') as f:
                f.write('v3')

            archive = _promote_v3_files({v3: orig}, genome_dir=tmp, genome_id='TEST')

            with tarfile.open(archive, 'r:gz') as tar:
                self.assertIn('2_cds/genome.gbk', tar.getnames())


class TestNormalizeLtDigits(TestCase):
    """normalize() must produce 6-digit locus tags to match the v3 spec."""

    def test_normalize_produces_6_digit_locus_tags(self):
        with (tempfile.NamedTemporaryFile(suffix='.gbk', delete=False) as gbk_f,
              tempfile.NamedTemporaryFile(suffix='.gbk', delete=False) as out_f):
            gbk_path, out_path = gbk_f.name, out_f.name

        try:
            _write_gbk(gbk_path, [('old_contig', ['OLD_00001', 'OLD_00002', 'OLD_00003'])])
            contig_map, lt_map = GenBankFile(gbk_path).normalize(out=out_path, genome_id=GENOME_ID)
            new_lts = list(lt_map.values())
            self.assertEqual(new_lts[0], f'{GENOME_ID}_000001')
            self.assertEqual(new_lts[1], f'{GENOME_ID}_000002')
            # v3 checker must accept the output
            genome_json = {'identifier': GENOME_ID, 'cds_tool_gbk_file': os.path.basename(out_path)}
            genome_dir = os.path.dirname(out_path)
            with open(os.path.join(genome_dir, 'genome.json'), 'w') as f:
                json.dump(genome_json, f)
            result = check_genome_v3(genome_dir, GENOME_ID)
            self.assertTrue(result.is_v3, f'post-normalize check failed: {result.issues}')
        finally:
            for p in (gbk_path, out_path):
                if os.path.exists(p):
                    os.unlink(p)
            genome_json_path = os.path.join(os.path.dirname(out_path), 'genome.json')
            if os.path.exists(genome_json_path):
                os.unlink(genome_json_path)


class TestFullUpgradeFlow(TestCase):
    """End-to-end test of the v2→v3 upgrade on a single genome dir (bypassing folder_structure)."""

    def _setup_v2_genome(self, tmp: str) -> str:
        """Create a minimal v2 genome folder and return its path."""
        genome_dir = os.path.join(tmp, 'organisms', GENOME_ID, 'genomes', GENOME_ID)
        os.makedirs(genome_dir)

        gbk_path = os.path.join(genome_dir, f'{GENOME_ID}.gbk')
        fna_path = os.path.join(genome_dir, f'{GENOME_ID}.fna')
        ann_path = os.path.join(genome_dir, f'{GENOME_ID}.KG')

        _write_gbk(gbk_path, [('old_contig_1', ['OLD_00001', 'OLD_00002']),
                               ('old_contig_2', ['OLD_00003'])])
        _write_fna(fna_path, ['old_contig_1', 'old_contig_2'])
        _write_annotation(ann_path, [('OLD_00001', 'K00001'), ('OLD_00002', 'K00002'),
                                     ('OLD_00003', 'K00003')])

        genome_json = {
            'identifier': GENOME_ID,
            'cds_tool_gbk_file': f'{GENOME_ID}.gbk',
            'assembly_fasta_file': f'{GENOME_ID}.fna',
            'custom_annotations': [{'file': f'{GENOME_ID}.KG', 'type': 'KG'}],
        }
        with open(os.path.join(genome_dir, 'genome.json'), 'w') as f:
            json.dump(genome_json, f)

        return genome_dir

    def test_pre_check_not_v3(self):
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._setup_v2_genome(tmp)
            result = check_genome_v3(genome_dir, GENOME_ID)
            self.assertFalse(result.is_v3)

    def test_normalize_then_check_v3(self):
        """Normalize the GBK, update FNA headers, update annotations → post-check passes."""
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._setup_v2_genome(tmp)
            gbk_path = os.path.join(genome_dir, f'{GENOME_ID}.gbk')
            fna_path = os.path.join(genome_dir, f'{GENOME_ID}.fna')
            ann_path = os.path.join(genome_dir, f'{GENOME_ID}.KG')

            # --- Generate .v3 files ---
            from arx_tools.update_folder_structure import _apply_lt_map_to_file, _apply_contig_map_to_fna, _promote_v3_files

            gbk_v3 = gbk_path + '.v3'
            contig_map, lt_map = GenBankFile(gbk_path).normalize(out=gbk_v3, genome_id=GENOME_ID)

            # Assembly FNA
            asm_v3 = fna_path + '.v3'
            _apply_contig_map_to_fna(fna_path, asm_v3, contig_map)

            # Annotation
            ann_v3 = ann_path + '.v3'
            _apply_lt_map_to_file(ann_path, ann_v3, lt_map)

            v3_to_orig = {
                gbk_v3: gbk_path,
                asm_v3: fna_path,
                ann_v3: ann_path,
            }
            _promote_v3_files(v3_to_orig, genome_dir=genome_dir, genome_id=GENOME_ID)

            # --- Post-check ---
            result = check_genome_v3(genome_dir, GENOME_ID, deep=True)
            self.assertTrue(result.is_v3, f'Post-upgrade check failed: {result.issues}')

            # --- Originals archived ---
            archive = os.path.join(genome_dir, f'{GENOME_ID}_v2_backup.tar.gz')
            self.assertTrue(os.path.exists(archive))
            with tarfile.open(archive, 'r:gz') as tar:
                archived = tar.getnames()
            self.assertTrue(any('gbk' in n for n in archived))
            self.assertTrue(any('fna' in n for n in archived))
            self.assertTrue(any('.KG' in n for n in archived))

    def test_failed_normalization_cleans_up_v3_files(self):
        """
        If GBK normalisation raises, from_2_to_3 must delete any .v3 files it created
        and leave originals untouched.
        """
        from unittest.mock import patch
        from arx_tools import update_folder_structure as ufs

        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._setup_v2_genome(tmp)

            # Write a version.json so ask() / set_folder_structure_version() work
            version_json = os.path.join(tmp, 'version.json')
            with open(version_json, 'w') as f:
                json.dump({'folder_structure_version': 2}, f)

            gbk_path = os.path.join(genome_dir, f'{GENOME_ID}.gbk')
            original_gbk = open(gbk_path).read()

            # Patch normalize() to raise after creating the .v3 file, and bypass the
            # interactive prompt and folder-structure version bump.
            def bad_normalize(_, out, **__):
                # Create the output file (as the real normalize would start to) then blow up
                with open(out, 'w') as f:
                    f.write('partial')
                raise RuntimeError('simulated disk error')

            with (patch.object(GenBankFile, 'normalize', bad_normalize),
                  patch.object(ufs, 'ask'),
                  patch.object(ufs, 'set_folder_structure_version')):
                ufs.from_2_to_3(folder_structure_dir=tmp)

            # Original GBK untouched
            with open(gbk_path) as f:
                self.assertEqual(f.read(), original_gbk)
            # .v3 file cleaned up
            self.assertFalse(os.path.exists(gbk_path + '.v3'))
            # No archive created
            self.assertFalse(os.path.exists(
                os.path.join(genome_dir, f'{GENOME_ID}_v2_backup.tar.gz')))
