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
    _apply_gene_tag_map_to_file,
    _apply_gene_tag_map_to_fasta,
    _apply_contig_map_to_fna,
    _apply_maps_to_gff,
    _promote_v3_files,
    _extend_gene_tag_map,
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
            _apply_gene_tag_map_to_file(src_path, dst_path, lt_map)
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
            _apply_gene_tag_map_to_file(src_path, dst_path, lt_map, is_eggnog=True)
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
            _apply_gene_tag_map_to_file(src_path, dst_path, lt_map)
            with open(dst_path) as f:
                content = f.read()
            self.assertIn('NEW_000001', content)
            self.assertIn('UNKNOWN_TAG', content)
        finally:
            os.unlink(src_path)
            if os.path.exists(dst_path):
                os.unlink(dst_path)

    def test_identity_map_does_not_warn(self):
        """
        When locus tags are already v3 but contig IDs still need renaming,
        normalize() produces an identity gene_tag_map.  The annotation file's tags
        are found in the map (even though nothing changes), so matched > 0 and
        no spurious 'no locus tags matched' warning is triggered.
        """
        lt_map = {'OLD_000001': 'OLD_000001'}  # identity
        with tempfile.NamedTemporaryFile(mode='w', suffix='.KG', delete=False) as src:
            src.write('OLD_000001\tK00001\nOLD_000002\tK00002\n')
            src_path = src.name
        dst_path = src_path + '.out'
        try:
            matched, total = _apply_gene_tag_map_to_file(src_path, dst_path, lt_map)
            self.assertEqual(total, 2)
            self.assertEqual(matched, 1)  # OLD_000001 found; OLD_000002 not in map
        finally:
            os.unlink(src_path)
            if os.path.exists(dst_path):
                os.unlink(dst_path)

    def test_source_file_unchanged(self):
        """src is never modified: only dst is written."""
        lt_map = {'OLD_000001': 'NEW_000001'}
        original = 'OLD_000001\tK00001\n'
        with tempfile.NamedTemporaryFile(mode='w', suffix='.KG', delete=False) as src:
            src.write(original)
            src_path = src.name
        dst_path = src_path + '.out'
        try:
            _apply_gene_tag_map_to_file(src_path, dst_path, lt_map)
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
            from arx_tools.update_folder_structure import _apply_gene_tag_map_to_file as _apply_gene_tag_map_to_file, _apply_contig_map_to_fna, _promote_v3_files

            gbk_v3 = gbk_path + '.v3'
            contig_map, lt_map = GenBankFile(gbk_path).normalize(out=gbk_v3, genome_id=GENOME_ID)

            # Assembly FNA
            asm_v3 = fna_path + '.v3'
            _apply_contig_map_to_fna(fna_path, asm_v3, contig_map)

            # Annotation
            ann_v3 = ann_path + '.v3'
            _apply_gene_tag_map_to_file(ann_path, ann_v3, lt_map)

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

    def test_failed_fna_update_cleans_up_all_v3_files(self):
        """
        If assembly FNA update raises after GBK .v3 was already written,
        both .v3 files are cleaned up and originals are untouched.
        """
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
        """
        If annotation update raises after GBK .v3 was already written,
        all .v3 files are cleaned up and originals are untouched.
        """
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


class TestCheckGenomeV3(TestCase):
    def _setup_genome_dir(self, tmp: str, contig_ids: list[str], locus_tags_per_contig: list[list[str]],
                          fna_contig_ids: list[str] = None, custom_annotations: list = None) -> str:
        genome_dir = os.path.join(tmp, GENOME_ID)
        os.makedirs(genome_dir)
        gbk_path = os.path.join(genome_dir, f'{GENOME_ID}.gbk')
        fna_path = os.path.join(genome_dir, f'{GENOME_ID}.fna')
        _write_gbk(gbk_path, list(zip(contig_ids, locus_tags_per_contig)))
        _write_fna(fna_path, fna_contig_ids if fna_contig_ids is not None else contig_ids)
        genome_json = {
            'identifier': GENOME_ID,
            'cds_tool_gbk_file': f'{GENOME_ID}.gbk',
            'assembly_fasta_file': f'{GENOME_ID}.fna',
            'custom_annotations': custom_annotations or [],
        }
        with open(os.path.join(genome_dir, 'genome.json'), 'w') as f:
            json.dump(genome_json, f)
        return genome_dir

    def _v3_genome_dir(self, tmp: str) -> str:
        return self._setup_genome_dir(
            tmp,
            contig_ids=[f'{GENOME_ID}_scf1', f'{GENOME_ID}_scf2'],
            locus_tags_per_contig=[[f'{GENOME_ID}_000001', f'{GENOME_ID}_000002'], [f'{GENOME_ID}_000003']],
        )

    def _v2_genome_dir(self, tmp: str) -> str:
        return self._setup_genome_dir(
            tmp,
            contig_ids=['old_contig_1', 'old_contig_2'],
            locus_tags_per_contig=[['OLD_00001', 'OLD_00002'], ['OLD_00003']],
        )

    def test_v3_genome_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._v3_genome_dir(tmp)
            result = check_genome_v3(genome_dir, GENOME_ID)
            self.assertTrue(result.is_v3)
            self.assertFalse(result.has_pending_v3_files)
            self.assertEqual(result.issues, [])

    def test_bad_gbk_contig_ids_detected(self):
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._v2_genome_dir(tmp)
            result = check_genome_v3(genome_dir, GENOME_ID)
            self.assertFalse(result.is_v3)
            self.assertTrue(any('contig' in i.lower() for i in result.issues))

    def test_bad_gbk_locus_tags_detected(self):
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._v2_genome_dir(tmp)
            result = check_genome_v3(genome_dir, GENOME_ID)
            self.assertFalse(result.is_v3)
            self.assertTrue(any('locus_tag' in i for i in result.issues))

    def test_bad_fna_headers_detected(self):
        """FNA with non-v3 headers fails even when GBK is already v3."""
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._setup_genome_dir(
                tmp,
                contig_ids=[f'{GENOME_ID}_scf1'],
                locus_tags_per_contig=[[f'{GENOME_ID}_000001']],
                fna_contig_ids=['old_contig_1'],
            )
            result = check_genome_v3(genome_dir, GENOME_ID)
            self.assertFalse(result.is_v3)
            self.assertTrue(any('FNA' in i for i in result.issues))

    def test_pending_v3_files_detected_shallow(self):
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._v3_genome_dir(tmp)
            gbk_path = os.path.join(genome_dir, f'{GENOME_ID}.gbk')
            with open(gbk_path + '.v3', 'w') as f:
                f.write('partial')
            result = check_genome_v3(genome_dir, GENOME_ID)
            self.assertTrue(result.has_pending_v3_files)
            self.assertIn(gbk_path + '.v3', result.pending_files)

    def test_pending_v3_annotation_detected_in_deep_mode(self):
        with tempfile.TemporaryDirectory() as tmp:
            ann_file = f'{GENOME_ID}.KG'
            genome_dir = self._setup_genome_dir(
                tmp,
                contig_ids=[f'{GENOME_ID}_scf1'],
                locus_tags_per_contig=[[f'{GENOME_ID}_000001']],
                custom_annotations=[{'file': ann_file, 'type': 'KG'}],
            )
            ann_path = os.path.join(genome_dir, ann_file)
            _write_annotation(ann_path, [(f'{GENOME_ID}_000001', 'K00001')])
            with open(ann_path + '.v3', 'w') as f:
                f.write('partial')
            shallow = check_genome_v3(genome_dir, GENOME_ID, deep=False)
            self.assertFalse(shallow.has_pending_v3_files)
            deep = check_genome_v3(genome_dir, GENOME_ID, deep=True)
            self.assertTrue(deep.has_pending_v3_files)
            self.assertIn(ann_path + '.v3', deep.pending_files)

    def test_deep_check_detects_annotation_issues(self):
        """Shallow check passes for a v3 GBK/FNA; deep check catches non-v3 annotation tags."""
        with tempfile.TemporaryDirectory() as tmp:
            ann_file = f'{GENOME_ID}.KG'
            genome_dir = self._setup_genome_dir(
                tmp,
                contig_ids=[f'{GENOME_ID}_scf1'],
                locus_tags_per_contig=[[f'{GENOME_ID}_000001']],
                custom_annotations=[{'file': ann_file, 'type': 'KG'}],
            )
            ann_path = os.path.join(genome_dir, ann_file)
            _write_annotation(ann_path, [('OLD_00001', 'K00001')])
            shallow = check_genome_v3(genome_dir, GENOME_ID, deep=False)
            self.assertTrue(shallow.is_v3)
            deep = check_genome_v3(genome_dir, GENOME_ID, deep=True)
            self.assertFalse(deep.is_v3)
            self.assertTrue(any(ann_file in i for i in deep.issues))

    def test_missing_genome_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = os.path.join(tmp, GENOME_ID)
            os.makedirs(genome_dir)
            result = check_genome_v3(genome_dir, GENOME_ID)
            self.assertFalse(result.is_v3)
            self.assertIn('genome.json not found', result.issues)

    def test_missing_gbk_on_disk(self):
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._v3_genome_dir(tmp)
            os.remove(os.path.join(genome_dir, f'{GENOME_ID}.gbk'))
            result = check_genome_v3(genome_dir, GENOME_ID)
            self.assertFalse(result.is_v3)
            self.assertTrue(any('GBK not found' in i for i in result.issues))

    def test_summary_v3_ok(self):
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._v3_genome_dir(tmp)
            result = check_genome_v3(genome_dir, GENOME_ID)
            self.assertIn('v3 OK', result.summary(GENOME_ID))
            self.assertIn(GENOME_ID, result.summary(GENOME_ID))

    def test_summary_not_v3(self):
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._v2_genome_dir(tmp)
            result = check_genome_v3(genome_dir, GENOME_ID)
            summary = result.summary(GENOME_ID)
            self.assertIn('NOT v3', summary)
            self.assertIn(GENOME_ID, summary)

    def test_summary_partial_upgrade(self):
        with tempfile.TemporaryDirectory() as tmp:
            genome_dir = self._v3_genome_dir(tmp)
            gbk_path = os.path.join(genome_dir, f'{GENOME_ID}.gbk')
            with open(gbk_path + '.v3', 'w') as f:
                f.write('partial')
            result = check_genome_v3(genome_dir, GENOME_ID)
            self.assertIn('PARTIAL', result.summary(GENOME_ID))


class TestFrom2To3SkipPaths(TestCase):
    def _setup_fs(self, tmp: str, gbk_contig_ids: list[str] = None,
                  gbk_locus_tags: list[list[str]] = None) -> str:
        """Create a minimal folder structure with one genome and version.json."""
        os.makedirs(os.path.join(tmp, 'organisms', GENOME_ID, 'genomes', GENOME_ID))
        genome_dir = os.path.join(tmp, 'organisms', GENOME_ID, 'genomes', GENOME_ID)
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
        """If a .v3 file already exists the genome is warned and skipped; originals are untouched."""
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


class TestApplyMapsToGff(TestCase):
    def _write(self, path, lines):
        with open(path, 'w') as f:
            f.writelines(lines)

    def _roundtrip(self, lines, contig_map, gene_tag_map):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as src:
            src.writelines(lines)
            src_path = src.name
        dst_path = src_path + '.out'
        try:
            changed = _apply_maps_to_gff(src_path, dst_path, contig_map, gene_tag_map)
            with open(dst_path) as f:
                content = f.read()
            return content, changed
        finally:
            os.unlink(src_path)
            if os.path.exists(dst_path):
                os.unlink(dst_path)

    def test_data_line_contig_and_locus_tag_renamed(self):
        line = 'old_scf1\t.\tCDS\t1\t9\t.\t+\t0\tID=OLD_000001;locus_tag=OLD_000001\n'
        content, changed = self._roundtrip(
            [line],
            contig_map={'old_scf1': 'NEW_scf1'},
            gene_tag_map={'OLD_000001': 'NEW_000001'},
        )
        self.assertIn('NEW_scf1\t', content)
        self.assertIn('NEW_000001', content)
        self.assertNotIn('old_scf1', content)
        self.assertNotIn('OLD_000001', content)
        self.assertEqual(changed, 1)

    def test_sequence_region_pragma_renamed(self):
        lines = ['##gff-version 3\n', '##sequence-region old_scf1 1 1000\n']
        content, _ = self._roundtrip(lines, contig_map={'old_scf1': 'NEW_scf1'}, gene_tag_map={})
        self.assertIn('##sequence-region NEW_scf1 1 1000', content)
        self.assertNotIn('old_scf1', content)

    def test_embedded_fasta_contig_renamed(self):
        lines = ['##gff-version 3\n', '##FASTA\n', '>old_scf1 desc\n', 'ATCG\n']
        content, _ = self._roundtrip(lines, contig_map={'old_scf1': 'NEW_scf1'}, gene_tag_map={})
        self.assertIn('>NEW_scf1 desc', content)
        self.assertNotIn('>old_scf1', content)
        self.assertIn('ATCG', content)

    def test_unmapped_values_preserved(self):
        line = 'old_scf1\t.\tCDS\t1\t9\t.\t+\t0\tID=OLD_000001\n'
        content, changed = self._roundtrip(
            [line],
            contig_map={'other_scf': 'X'},
            gene_tag_map={'OTHER_000001': 'Y'},
        )
        self.assertIn('old_scf1', content)
        self.assertIn('OLD_000001', content)
        self.assertEqual(changed, 0)

    def test_short_lines_passed_through_unchanged(self):
        """Lines with fewer than 9 tab-separated columns are not modified."""
        line = 'old_scf1\t.\tgene\n'
        content, _ = self._roundtrip([line], contig_map={'old_scf1': 'NEW_scf1'}, gene_tag_map={})
        self.assertEqual(content, line)

    def test_source_file_unchanged(self):
        original = 'old_scf1\t.\tCDS\t1\t9\t.\t+\t0\tID=foo\n'
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as src:
            src.write(original)
            src_path = src.name
        dst_path = src_path + '.out'
        try:
            _apply_maps_to_gff(src_path, dst_path, {'old_scf1': 'NEW_scf1'}, {})
            with open(src_path) as f:
                self.assertEqual(f.read(), original)
        finally:
            os.unlink(src_path)
            if os.path.exists(dst_path):
                os.unlink(dst_path)


class TestApplyGeneTagMapToFasta(TestCase):
    def _roundtrip(self, content, gene_tag_map):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.faa', delete=False) as src:
            src.write(content)
            src_path = src.name
        dst_path = src_path + '.out'
        try:
            count = _apply_gene_tag_map_to_fasta(src_path, dst_path, gene_tag_map)
            with open(dst_path) as f:
                result = f.read()
            return result, count
        finally:
            os.unlink(src_path)
            if os.path.exists(dst_path):
                os.unlink(dst_path)

    def test_basic_rename(self):
        result, count = self._roundtrip(
            '>OLD_000001 some product\nMPKL\n',
            {'OLD_000001': 'NEW_000001'},
        )
        self.assertEqual(count, 1)
        self.assertIn('>NEW_000001 some product', result)
        self.assertNotIn('>OLD_', result)

    def test_description_preserved(self):
        result, _ = self._roundtrip(
            '>OLD_000001 hypothetical protein [Organism X]\nMPKL\n',
            {'OLD_000001': 'NEW_000001'},
        )
        self.assertIn('hypothetical protein [Organism X]', result)

    def test_unmapped_tag_unchanged(self):
        result, count = self._roundtrip('>OLD_000001\nMPKL\n', {'OTHER': 'NEW'})
        self.assertEqual(count, 0)
        self.assertIn('>OLD_000001', result)

    def test_sequence_lines_not_modified(self):
        """Sequence lines that look like a locus tag are never touched."""
        result, _ = self._roundtrip(
            '>OLD_000001\nOLD_000001\n',
            {'OLD_000001': 'NEW_000001'},
        )
        lines = result.splitlines()
        self.assertEqual(lines[1], 'OLD_000001')

    def test_source_file_unchanged(self):
        original = '>OLD_000001\nMPKL\n'
        with tempfile.NamedTemporaryFile(mode='w', suffix='.faa', delete=False) as src:
            src.write(original)
            src_path = src.name
        dst_path = src_path + '.out'
        try:
            _apply_gene_tag_map_to_fasta(src_path, dst_path, {'OLD_000001': 'NEW_000001'})
            with open(src_path) as f:
                self.assertEqual(f.read(), original)
        finally:
            os.unlink(src_path)
            if os.path.exists(dst_path):
                os.unlink(dst_path)


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


class TestExtendGeneTagMap(TestCase):
    def test_adds_five_digit_variant(self):
        """5-digit variant of each v3 tag is added as a key mapping to the same v3 tag."""
        extended = _extend_gene_tag_map({'NCBI_RS001': f'{GENOME_ID}_000001'})
        self.assertEqual(extended['NCBI_RS001'], f'{GENOME_ID}_000001')
        self.assertEqual(extended.get(f'{GENOME_ID}_00001'), f'{GENOME_ID}_000001')

    def test_no_extra_entry_when_five_digit_already_present(self):
        """If the 5-digit key already exists (Prokka-style GBK), no extra entry is added."""
        gene_map = {f'{GENOME_ID}_00001': f'{GENOME_ID}_000001'}
        extended = _extend_gene_tag_map(gene_map)
        self.assertEqual(len(extended), len(gene_map))

    def test_number_above_99999_not_added(self):
        """Numbers > 99999 already occupy 6 digits; no 5-digit variant is generated."""
        extended = _extend_gene_tag_map({'OLD': f'{GENOME_ID}_100001'})
        self.assertNotIn(f'{GENOME_ID}_100001', set(extended) - {'OLD'})

    def test_empty_map_returns_empty(self):
        self.assertEqual(_extend_gene_tag_map({}), {})

    def test_complex_genome_id_prefix(self):
        """Genome IDs containing dots and hyphens (e.g. NCBI accession style) are handled."""
        genome_id = 'MOD1-EC5552_GCF_002228865.1_ASM222886v1'
        extended = _extend_gene_tag_map({'BEG76_RS26710': f'{genome_id}_000001'})
        self.assertEqual(extended.get(f'{genome_id}_00001'), f'{genome_id}_000001')


class TestArxAssignedLocusTagsInAnnotations(TestCase):
    """
    When the GBK holds external locus tags (e.g. NCBI RefSeq IDs) but annotation
    files and derived FASTA files were generated using arx-assigned 5-digit locus
    tags, the extended map must rename those tags to v3 6-digit format.
    """

    def _make_file(self, content: str, suffix: str) -> str:
        f = tempfile.NamedTemporaryFile(mode='w', suffix=suffix, delete=False)
        f.write(content)
        f.close()
        return f.name

    def _roundtrip(self, src_path: str, lt_map: dict, is_eggnog: bool = False) -> str:
        dst_path = src_path + '.out'
        try:
            _apply_gene_tag_map_to_file(src_path, dst_path, lt_map, is_eggnog=is_eggnog)
            with open(dst_path) as f:
                return f.read()
        finally:
            if os.path.exists(dst_path):
                os.unlink(dst_path)

    def test_arx_five_digit_in_regular_annotation(self):
        """AR-type (and other regular) annotation files with 5-digit arx locus tags are renamed."""
        lt_map = _extend_gene_tag_map({'NCBI_TAG': f'{GENOME_ID}_000001'})
        src = self._make_file(f'{GENOME_ID}_00001\tAR:blaNDM-1\n', '.AR')
        try:
            content = self._roundtrip(src, lt_map)
            self.assertIn(f'{GENOME_ID}_000001', content)
            self.assertNotIn(f'{GENOME_ID}_00001', content)
        finally:
            os.unlink(src)

    def test_arx_five_digit_in_eggnog_with_pipe_prefix(self):
        """Old-style eggnog files with gnl|extdb| prefix and 5-digit arx locus tags are renamed."""
        lt_map = _extend_gene_tag_map({'NCBI_TAG': f'{GENOME_ID}_000001'})
        src = self._make_file(f'#query_name\theader\ngnl|extdb|{GENOME_ID}_00001\tdata\n', '.annotations')
        try:
            content = self._roundtrip(src, lt_map, is_eggnog=True)
            self.assertIn(f'gnl|extdb|{GENOME_ID}_000001', content)
            self.assertNotIn(f'{GENOME_ID}_00001', content)
        finally:
            os.unlink(src)

    def test_arx_five_digit_in_eggnog_without_pipe(self):
        """eggnog-2.1.2 files with bare 5-digit arx locus tags (no pipe) are renamed."""
        lt_map = _extend_gene_tag_map({'NCBI_TAG': f'{GENOME_ID}_000001'})
        src = self._make_file(f'#query\theader\n{GENOME_ID}_00001\tdata\n', '.annotations')
        try:
            content = self._roundtrip(src, lt_map, is_eggnog=True)
            self.assertIn(f'{GENOME_ID}_000001', content)
            self.assertNotIn(f'{GENOME_ID}_00001\t', content)
        finally:
            os.unlink(src)

    def test_arx_five_digit_in_fasta(self):
        """FAA/FFN files with 5-digit arx locus tags in headers are renamed."""
        lt_map = _extend_gene_tag_map({'NCBI_TAG': f'{GENOME_ID}_000001'})
        with tempfile.NamedTemporaryFile(mode='w', suffix='.faa', delete=False) as f:
            f.write(f'>{GENOME_ID}_00001 some product\nMPKL\n')
            src_path = f.name
        dst_path = src_path + '.out'
        try:
            count = _apply_gene_tag_map_to_fasta(src_path, dst_path, lt_map)
            self.assertEqual(count, 1)
            with open(dst_path) as f:
                content = f.read()
            self.assertIn(f'>{GENOME_ID}_000001', content)
            self.assertNotIn(f'>{GENOME_ID}_00001', content)
        finally:
            os.unlink(src_path)
            if os.path.exists(dst_path):
                os.unlink(dst_path)

    def test_no_false_match_when_prokka_style_gbk(self):
        """When GBK already uses arx-style locus tags, direct lookup still works."""
        gene_map = {f'{GENOME_ID}_00001': f'{GENOME_ID}_000001'}
        lt_map = _extend_gene_tag_map(gene_map)
        src = self._make_file(f'{GENOME_ID}_00001\tK00001\n', '.KG')
        try:
            content = self._roundtrip(src, lt_map)
            self.assertIn(f'{GENOME_ID}_000001', content)
        finally:
            os.unlink(src)
