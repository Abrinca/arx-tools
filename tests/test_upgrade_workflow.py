"""
Workflow-level tests for the v2→v3 folder-structure upgrade.

Covers:
  _promote_v3_files
  GenBankFile.normalize() → 6-digit locus tags
  from_2_to_3() end-to-end: success, skip, failure-rollback, version bump
  --create_only / --promote two-step workflow
"""
import json
import os
import tarfile
import tempfile
from unittest import TestCase

from Bio import SeqIO

from arx_tools.update_folder_structure import (
    _apply_contig_map_to_fna,
    _apply_gene_tag_map_to_file,
    _promote_v3_files,
)
from arx_tools.check_v3 import check_genome_v3
from arx_tools.rename_genbank import GenBankFile

from tests.helpers import GENOME_ID, _write_gbk, _write_fna, _write_annotation


# ── _promote_v3_files ─────────────────────────────────────────────────────────

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


# ── normalize() 6-digit locus tags ───────────────────────────────────────────

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


# ── Full upgrade flow ─────────────────────────────────────────────────────────

class TestFullUpgradeFlow(TestCase):
    """End-to-end test of the v2→v3 upgrade on a single genome dir."""

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

            gbk_v3 = gbk_path + '.v3'
            contig_map, lt_map = GenBankFile(gbk_path).normalize(out=gbk_v3, genome_id=GENOME_ID)

            asm_v3 = fna_path + '.v3'
            _apply_contig_map_to_fna(fna_path, asm_v3, contig_map)

            ann_v3 = ann_path + '.v3'
            _apply_gene_tag_map_to_file(ann_path, ann_v3, lt_map)

            v3_to_orig = {gbk_v3: gbk_path, asm_v3: fna_path, ann_v3: ann_path}
            _promote_v3_files(v3_to_orig, genome_dir=genome_dir, genome_id=GENOME_ID)

            result = check_genome_v3(genome_dir, GENOME_ID, deep=True)
            self.assertTrue(result.is_v3, f'Post-upgrade check failed: {result.issues}')

            archive = os.path.join(genome_dir, f'{GENOME_ID}_v2_backup.tar.gz')
            self.assertTrue(os.path.exists(archive))
            with tarfile.open(archive, 'r:gz') as tar:
                archived = tar.getnames()
            self.assertTrue(any('gbk' in n for n in archived))
            self.assertTrue(any('fna' in n for n in archived))
            self.assertTrue(any('.KG' in n for n in archived))

    def test_failed_normalization_cleans_up_v3_files(self):
        """If GBK normalisation raises, from_2_to_3 must delete any .v3 files it created."""
        from unittest.mock import patch
        from arx_tools import update_folder_structure as ufs

        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._setup_v2_genome(tmp)

            with open(os.path.join(tmp, 'version.json'), 'w') as f:
                json.dump({'folder_structure_version': 2}, f)

            gbk_path = os.path.join(genome_dir, f'{GENOME_ID}.gbk')
            original_gbk = open(gbk_path).read()

            def bad_normalize(_, out, **__):
                with open(out, 'w') as f:
                    f.write('partial')
                raise RuntimeError('simulated disk error')

            with (patch.object(GenBankFile, 'normalize', bad_normalize),
                  patch.object(ufs, 'ask'),
                  patch.object(ufs, 'set_folder_structure_version')):
                ufs.from_2_to_3(folder_structure_dir=tmp)

            with open(gbk_path) as f:
                self.assertEqual(f.read(), original_gbk)
            self.assertFalse(os.path.exists(gbk_path + '.v3'))
            self.assertFalse(os.path.exists(
                os.path.join(genome_dir, f'{GENOME_ID}_v2_backup.tar.gz')))

    def test_failed_fna_update_cleans_up_all_v3_files(self):
        """If assembly FNA update raises after GBK .v3 was written, both .v3 files are removed."""
        from unittest.mock import patch
        from arx_tools import update_folder_structure as ufs

        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._setup_v2_genome(tmp)
            with open(os.path.join(tmp, 'version.json'), 'w') as f:
                json.dump({'folder_structure_version': 2}, f)

            gbk_path = os.path.join(genome_dir, f'{GENOME_ID}.gbk')
            fna_path = os.path.join(genome_dir, f'{GENOME_ID}.fna')
            original_fna = open(fna_path).read()

            def bad_fna(src, dst, contig_map):
                with open(dst, 'w') as f:
                    f.write('partial')
                raise RuntimeError('simulated fna error')

            with (patch.object(ufs, '_apply_contig_map_to_fna', bad_fna),
                  patch.object(ufs, 'ask'),
                  patch.object(ufs, 'set_folder_structure_version')):
                ufs.from_2_to_3(folder_structure_dir=tmp)

            with open(fna_path) as f:
                self.assertEqual(f.read(), original_fna)
            self.assertFalse(os.path.exists(gbk_path + '.v3'))
            self.assertFalse(os.path.exists(fna_path + '.v3'))
            self.assertFalse(os.path.exists(os.path.join(genome_dir, f'{GENOME_ID}_v2_backup.tar.gz')))

    def test_failed_annotation_update_cleans_up_all_v3_files(self):
        """If annotation update raises, all .v3 files are cleaned up."""
        from unittest.mock import patch
        from arx_tools import update_folder_structure as ufs

        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._setup_v2_genome(tmp)
            with open(os.path.join(tmp, 'version.json'), 'w') as f:
                json.dump({'folder_structure_version': 2}, f)

            gbk_path = os.path.join(genome_dir, f'{GENOME_ID}.gbk')
            ann_path = os.path.join(genome_dir, f'{GENOME_ID}.KG')
            original_ann = open(ann_path).read()

            def bad_annotation(src, dst, lt_map, **__):
                with open(dst, 'w') as f:
                    f.write('partial')
                raise RuntimeError('simulated annotation error')

            with (patch.object(ufs, '_apply_gene_tag_map_to_file', bad_annotation),
                  patch.object(ufs, 'ask'),
                  patch.object(ufs, 'set_folder_structure_version')):
                ufs.from_2_to_3(folder_structure_dir=tmp)

            with open(ann_path) as f:
                self.assertEqual(f.read(), original_ann)
            self.assertFalse(os.path.exists(gbk_path + '.v3'))
            self.assertFalse(os.path.exists(ann_path + '.v3'))
            self.assertFalse(os.path.exists(os.path.join(genome_dir, f'{GENOME_ID}_v2_backup.tar.gz')))


# ── from_2_to_3 skip paths ────────────────────────────────────────────────────

class TestFrom2To3SkipPaths(TestCase):
    def _setup_fs(self, tmp: str, gbk_contig_ids: list[str] = None,
                  gbk_locus_tags: list[list[str]] = None) -> str:
        """Create a minimal folder structure with one genome and version.json."""
        genome_dir = os.path.join(tmp, 'organisms', GENOME_ID, 'genomes', GENOME_ID)
        os.makedirs(genome_dir)
        gbk_path = os.path.join(genome_dir, f'{GENOME_ID}.gbk')
        fna_path = os.path.join(genome_dir, f'{GENOME_ID}.fna')
        contig_ids = gbk_contig_ids or ['old_contig_1']
        locus_tags = gbk_locus_tags or [['OLD_00001', 'OLD_00002']]
        _write_gbk(gbk_path, list(zip(contig_ids, locus_tags)))
        _write_fna(fna_path, contig_ids)
        genome_json = {
            'identifier': GENOME_ID,
            'cds_tool_gbk_file': f'{GENOME_ID}.gbk',
            'assembly_fasta_file': f'{GENOME_ID}.fna',
            'custom_annotations': [],
        }
        with open(os.path.join(genome_dir, 'genome.json'), 'w') as f:
            json.dump(genome_json, f)
        with open(os.path.join(tmp, 'version.json'), 'w') as f:
            json.dump({'folder_structure_version': 2}, f)
        return genome_dir

    def _run(self, tmp: str) -> None:
        from unittest.mock import patch
        from arx_tools import update_folder_structure as ufs
        with (patch.object(ufs, 'ask'),
              patch.object(ufs, 'set_folder_structure_version')):
            ufs.from_2_to_3(folder_structure_dir=tmp)

    def test_already_v3_is_skipped(self):
        """A genome whose GBK/FNA already use v3 IDs is skipped without creating any files."""
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._setup_fs(
                tmp,
                gbk_contig_ids=[f'{GENOME_ID}_scf1'],
                gbk_locus_tags=[[f'{GENOME_ID}_000001', f'{GENOME_ID}_000002']],
            )
            fna_path = os.path.join(genome_dir, f'{GENOME_ID}.fna')
            _write_fna(fna_path, [f'{GENOME_ID}_scf1'])
            self._run(tmp)
            for fname in os.listdir(genome_dir):
                self.assertFalse(fname.endswith('.v3'), f'Unexpected .v3 file: {fname}')
            self.assertFalse(os.path.exists(os.path.join(genome_dir, f'{GENOME_ID}_v2_backup.tar.gz')))

    def test_partial_upgrade_is_skipped(self):
        """If a .v3 file already exists the genome is warned and skipped; originals untouched."""
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._setup_fs(tmp)
            gbk_path = os.path.join(genome_dir, f'{GENOME_ID}.gbk')
            original_content = open(gbk_path).read()
            with open(gbk_path + '.v3', 'w') as f:
                f.write('leftover')
            self._run(tmp)
            with open(gbk_path) as f:
                self.assertEqual(f.read(), original_content)
            self.assertFalse(os.path.exists(os.path.join(genome_dir, f'{GENOME_ID}_v2_backup.tar.gz')))

    def test_missing_cds_tool_gbk_file_is_skipped(self):
        """Genome with no cds_tool_gbk_file key in genome.json is skipped."""
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._setup_fs(tmp)
            json_path = os.path.join(genome_dir, 'genome.json')
            with open(json_path) as f:
                gj = json.load(f)
            del gj['cds_tool_gbk_file']
            with open(json_path, 'w') as f:
                json.dump(gj, f)
            self._run(tmp)
            self.assertFalse(os.path.exists(os.path.join(genome_dir, f'{GENOME_ID}_v2_backup.tar.gz')))

    def test_gbk_file_not_found_on_disk_is_skipped(self):
        """Genome whose GBK path in genome.json doesn't exist on disk is skipped."""
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._setup_fs(tmp)
            os.remove(os.path.join(genome_dir, f'{GENOME_ID}.gbk'))
            self._run(tmp)
            self.assertFalse(os.path.exists(os.path.join(genome_dir, f'{GENOME_ID}_v2_backup.tar.gz')))


# ── Version-bump behaviour ────────────────────────────────────────────────────

class TestVersionBumpBehavior(TestCase):
    def _setup_fs(self, tmp: str) -> str:
        genome_dir = os.path.join(tmp, 'organisms', GENOME_ID, 'genomes', GENOME_ID)
        os.makedirs(genome_dir)
        gbk_path = os.path.join(genome_dir, f'{GENOME_ID}.gbk')
        fna_path = os.path.join(genome_dir, f'{GENOME_ID}.fna')
        _write_gbk(gbk_path, [('old_contig_1', ['OLD_00001', 'OLD_00002'])])
        _write_fna(fna_path, ['old_contig_1'])
        genome_json = {
            'identifier': GENOME_ID,
            'cds_tool_gbk_file': f'{GENOME_ID}.gbk',
            'assembly_fasta_file': f'{GENOME_ID}.fna',
            'custom_annotations': [],
        }
        with open(os.path.join(genome_dir, 'genome.json'), 'w') as f:
            json.dump(genome_json, f)
        with open(os.path.join(tmp, 'version.json'), 'w') as f:
            json.dump({'folder_structure_version': 2}, f)
        return genome_dir

    def test_version_not_bumped_when_genome_fails(self):
        """set_folder_structure_version must NOT be called if any genome fails."""
        from unittest.mock import patch, MagicMock
        from arx_tools import update_folder_structure as ufs

        with tempfile.TemporaryDirectory() as tmp:
            self._setup_fs(tmp)

            def bad_normalize(self_gbk, out, **__):
                with open(out, 'w') as f:
                    f.write('partial')
                raise RuntimeError('simulated error')

            mock_set_version = MagicMock()
            with (patch.object(GenBankFile, 'normalize', bad_normalize),
                  patch.object(ufs, 'ask'),
                  patch.object(ufs, 'set_folder_structure_version', mock_set_version)):
                ufs.from_2_to_3(folder_structure_dir=tmp)

            mock_set_version.assert_not_called()

    def test_version_bumped_on_full_success(self):
        """version.json is updated to 3 when all genomes migrate without errors."""
        from unittest.mock import patch
        from arx_tools import update_folder_structure as ufs

        with tempfile.TemporaryDirectory() as tmp:
            self._setup_fs(tmp)
            with patch.object(ufs, 'ask'):
                ufs.from_2_to_3(folder_structure_dir=tmp)

            with open(os.path.join(tmp, 'version.json')) as f:
                version_data = json.load(f)
            self.assertEqual(version_data['folder_structure_version'], 3)


# ── --create_only / --promote two-step workflow ───────────────────────────────

class TestCreateOnlyAndPromote(TestCase):
    """Tests for the --create_only / --promote two-step upgrade workflow."""

    def _setup_fs(self, tmp: str) -> str:
        """Minimal v2 folder structure with one genome; returns the genome dir."""
        genome_dir = os.path.join(tmp, 'organisms', GENOME_ID, 'genomes', GENOME_ID)
        os.makedirs(genome_dir)
        gbk_path = os.path.join(genome_dir, f'{GENOME_ID}.gbk')
        fna_path = os.path.join(genome_dir, f'{GENOME_ID}.fna')
        ann_path = os.path.join(genome_dir, f'{GENOME_ID}.KG')
        _write_gbk(gbk_path, [('old_contig_1', ['OLD_00001', 'OLD_00002'])])
        _write_fna(fna_path, ['old_contig_1'])
        _write_annotation(ann_path, [('OLD_00001', 'K00001'), ('OLD_00002', 'K00002')])
        genome_json = {
            'identifier': GENOME_ID,
            'cds_tool_gbk_file': f'{GENOME_ID}.gbk',
            'assembly_fasta_file': f'{GENOME_ID}.fna',
            'custom_annotations': [{'file': f'{GENOME_ID}.KG', 'type': 'KG'}],
        }
        with open(os.path.join(genome_dir, 'genome.json'), 'w') as f:
            json.dump(genome_json, f)
        with open(os.path.join(tmp, 'version.json'), 'w') as f:
            json.dump({'folder_structure_version': 2}, f)
        return genome_dir

    def _run(self, tmp: str, **kwargs) -> None:
        from unittest.mock import patch
        from arx_tools import update_folder_structure as ufs
        with (patch.object(ufs, 'ask'),
              patch.object(ufs, 'set_folder_structure_version')):
            ufs.from_2_to_3(folder_structure_dir=tmp, **kwargs)

    # ── --create_only ─────────────────────────────────────────────────────────

    def test_create_only_creates_v3_files(self):
        """--create_only writes .v3 files for GBK, assembly FNA, and annotations."""
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._setup_fs(tmp)
            self._run(tmp, create_only=True)
            self.assertTrue(os.path.exists(os.path.join(genome_dir, f'{GENOME_ID}.gbk.v3')))
            self.assertTrue(os.path.exists(os.path.join(genome_dir, f'{GENOME_ID}.fna.v3')))
            self.assertTrue(os.path.exists(os.path.join(genome_dir, f'{GENOME_ID}.KG.v3')))

    def test_create_only_leaves_originals_untouched(self):
        """--create_only must not modify the original GBK or FNA."""
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._setup_fs(tmp)
            gbk_path = os.path.join(genome_dir, f'{GENOME_ID}.gbk')
            fna_path = os.path.join(genome_dir, f'{GENOME_ID}.fna')
            original_gbk = open(gbk_path).read()
            original_fna = open(fna_path).read()
            self._run(tmp, create_only=True)
            self.assertEqual(open(gbk_path).read(), original_gbk)
            self.assertEqual(open(fna_path).read(), original_fna)

    def test_create_only_does_not_create_archive(self):
        """--create_only must not create a backup archive."""
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._setup_fs(tmp)
            self._run(tmp, create_only=True)
            self.assertFalse(os.path.exists(
                os.path.join(genome_dir, f'{GENOME_ID}_v2_backup.tar.gz')))

    def test_create_only_does_not_bump_version(self):
        """--create_only must never call set_folder_structure_version."""
        from unittest.mock import patch, MagicMock
        from arx_tools import update_folder_structure as ufs
        with tempfile.TemporaryDirectory() as tmp:
            self._setup_fs(tmp)
            mock_set_version = MagicMock()
            with (patch.object(ufs, 'ask'),
                  patch.object(ufs, 'set_folder_structure_version', mock_set_version)):
                ufs.from_2_to_3(folder_structure_dir=tmp, create_only=True)
            mock_set_version.assert_not_called()

    def test_create_only_skips_genome_with_existing_pending_files(self):
        """If .v3 files already exist, --create_only warns and skips (no overwrite)."""
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._setup_fs(tmp)
            gbk_v3 = os.path.join(genome_dir, f'{GENOME_ID}.gbk.v3')
            with open(gbk_v3, 'w') as f:
                f.write('leftover')
            self._run(tmp, create_only=True)
            with open(gbk_v3) as f:
                self.assertEqual(f.read(), 'leftover')
            self.assertFalse(os.path.exists(
                os.path.join(genome_dir, f'{GENOME_ID}_v2_backup.tar.gz')))

    # ── --promote ─────────────────────────────────────────────────────────────

    def test_promote_promotes_pending_files(self):
        """--promote archives originals and moves .v3 files into place."""
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._setup_fs(tmp)
            gbk_path = os.path.join(genome_dir, f'{GENOME_ID}.gbk')
            with open(gbk_path + '.v3', 'w') as f:
                f.write('v3 content')
            self._run(tmp, promote=True)
            self.assertFalse(os.path.exists(gbk_path + '.v3'))
            with open(gbk_path) as f:
                self.assertEqual(f.read(), 'v3 content')
            self.assertTrue(os.path.exists(
                os.path.join(genome_dir, f'{GENOME_ID}_v2_backup.tar.gz')))

    def test_promote_skips_genome_without_pending_files(self):
        """--promote silently skips a genome that has no .v3 files."""
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._setup_fs(tmp)
            self._run(tmp, promote=True)
            self.assertFalse(os.path.exists(
                os.path.join(genome_dir, f'{GENOME_ID}_v2_backup.tar.gz')))

    def test_promote_does_not_bump_version_when_genomes_not_ready(self):
        """Version is not bumped when --promote finds genomes with no pending .v3 files."""
        from unittest.mock import patch, MagicMock
        from arx_tools import update_folder_structure as ufs
        with tempfile.TemporaryDirectory() as tmp:
            self._setup_fs(tmp)
            mock_set_version = MagicMock()
            with (patch.object(ufs, 'ask'),
                  patch.object(ufs, 'set_folder_structure_version', mock_set_version)):
                ufs.from_2_to_3(folder_structure_dir=tmp, promote=True)
            mock_set_version.assert_not_called()

    def test_promote_bumps_version_when_all_promoted(self):
        """Version is bumped to 3 after --promote succeeds for every genome."""
        from unittest.mock import patch
        from arx_tools import update_folder_structure as ufs
        with tempfile.TemporaryDirectory() as tmp:
            self._setup_fs(tmp)
            with patch.object(ufs, 'ask'):
                ufs.from_2_to_3(folder_structure_dir=tmp, create_only=True)
            with patch.object(ufs, 'ask'):
                ufs.from_2_to_3(folder_structure_dir=tmp, promote=True)
            with open(os.path.join(tmp, 'version.json')) as f:
                self.assertEqual(json.load(f)['folder_structure_version'], 3)

    # ── Two-step workflow ─────────────────────────────────────────────────────

    def test_create_only_then_promote_produces_v3_genome(self):
        """Full two-step workflow: --create_only followed by --promote yields a v3 genome."""
        from unittest.mock import patch
        from arx_tools import update_folder_structure as ufs
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._setup_fs(tmp)

            with patch.object(ufs, 'ask'):
                ufs.from_2_to_3(folder_structure_dir=tmp, create_only=True)
            pre = check_genome_v3(genome_dir, GENOME_ID)
            self.assertFalse(pre.is_v3, 'originals should still be v2 after --create_only')
            self.assertTrue(pre.has_pending_v3_files)

            with patch.object(ufs, 'ask'):
                ufs.from_2_to_3(folder_structure_dir=tmp, promote=True)
            post = check_genome_v3(genome_dir, GENOME_ID)
            self.assertTrue(post.is_v3, f'post-promote check failed: {post.issues}')
            self.assertFalse(post.has_pending_v3_files)

    def test_normal_run_does_not_bump_version_when_pending_files_exist(self):
        """A normal run that skips a genome due to leftover .v3 files must not bump the version."""
        from unittest.mock import patch, MagicMock
        from arx_tools import update_folder_structure as ufs
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._setup_fs(tmp)
            gbk_path = os.path.join(genome_dir, f'{GENOME_ID}.gbk')
            with open(gbk_path + '.v3', 'w') as f:
                f.write('leftover')
            mock_set_version = MagicMock()
            with (patch.object(ufs, 'ask'),
                  patch.object(ufs, 'set_folder_structure_version', mock_set_version)):
                ufs.from_2_to_3(folder_structure_dir=tmp)
            mock_set_version.assert_not_called()
