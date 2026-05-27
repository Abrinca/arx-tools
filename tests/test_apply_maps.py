"""
Unit tests for the individual map-application helpers in update_folder_structure,
and for import_genome normalization logic.

Covers:
  _apply_gene_tag_map_to_file
  _apply_gene_tag_map_to_fasta
  _apply_contig_map_to_fna
  _apply_maps_to_gff
  _extend_gene_tag_map
  _resolve_contig_id
  _lookup_gene_tag
  import_genome: gnl|C| prefix stripping before FNA/GBK contig ID comparison
"""
import os
import tempfile
from unittest import TestCase

from arx_tools.update_folder_structure import (
    _apply_gene_tag_map_to_file,
    _apply_gene_tag_map_to_fasta,
    _apply_contig_map_to_fna,
    _apply_maps_to_gff,
    _extend_gene_tag_map,
    _resolve_contig_id,
    _lookup_gene_tag,
)
from arx_tools.rename_fasta import FastaFile

from tests.helpers import GENOME_ID


# ── _apply_gene_tag_map_to_file ───────────────────────────────────────────────

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

    def test_identity_map_matched_count(self):
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


# ── _apply_contig_map_to_fna ──────────────────────────────────────────────────

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


class TestApplyContigMapToFnaGnlPrefix(TestCase):
    """_apply_contig_map_to_fna must rename gnl|C| prefixed headers (Prokka FNA)."""

    def test_gnl_prefix_renamed(self):
        contig_map = {'ALNJDMAK_1': 'GENOME_scf1', 'ALNJDMAK_2': 'GENOME_scf2'}
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fna', delete=False) as src:
            src.write('>gnl|C|ALNJDMAK_1 desc\nATCG\n>gnl|C|ALNJDMAK_2\nGGCC\n')
            src_path = src.name
        dst_path = src_path + '.out'
        try:
            count = _apply_contig_map_to_fna(src_path, dst_path, contig_map)
            self.assertEqual(count, 2)
            with open(dst_path) as f:
                content = f.read()
            self.assertIn('>GENOME_scf1 desc', content)
            self.assertIn('>GENOME_scf2\n', content)
            self.assertNotIn('gnl|C|', content)
        finally:
            os.unlink(src_path)
            if os.path.exists(dst_path):
                os.unlink(dst_path)

    def test_gnl_prefix_description_preserved(self):
        contig_map = {'ALNJDMAK_1': 'GENOME_scf1'}
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fna', delete=False) as src:
            src.write('>gnl|C|ALNJDMAK_1 length=40066 coverage=50x\nATCG\n')
            src_path = src.name
        dst_path = src_path + '.out'
        try:
            _apply_contig_map_to_fna(src_path, dst_path, contig_map)
            with open(dst_path) as f:
                content = f.read()
            self.assertIn('length=40066 coverage=50x', content)
        finally:
            os.unlink(src_path)
            if os.path.exists(dst_path):
                os.unlink(dst_path)


# ── _apply_maps_to_gff ────────────────────────────────────────────────────────

class TestApplyMapsToGff(TestCase):
    def _roundtrip(self, lines, contig_map, gene_tag_map):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as src:
            src.writelines(lines)
            src_path = src.name
        dst_path = src_path + '.out'
        try:
            seqid_changed, attr_renamed = _apply_maps_to_gff(src_path, dst_path, contig_map, gene_tag_map)
            with open(dst_path) as f:
                content = f.read()
            return content, seqid_changed, attr_renamed
        finally:
            os.unlink(src_path)
            if os.path.exists(dst_path):
                os.unlink(dst_path)

    def test_data_line_contig_and_locus_tag_renamed(self):
        line = 'old_scf1\t.\tCDS\t1\t9\t.\t+\t0\tID=OLD_000001;locus_tag=OLD_000001\n'
        content, seqid_changed, attr_renamed = self._roundtrip(
            [line],
            contig_map={'old_scf1': 'NEW_scf1'},
            gene_tag_map={'OLD_000001': 'NEW_000001'},
        )
        self.assertIn('NEW_scf1\t', content)
        self.assertIn('NEW_000001', content)
        self.assertNotIn('old_scf1', content)
        self.assertNotIn('OLD_000001', content)
        self.assertEqual(seqid_changed, 1)
        self.assertGreater(attr_renamed, 0)

    def test_sequence_region_pragma_renamed(self):
        lines = ['##gff-version 3\n', '##sequence-region old_scf1 1 1000\n']
        content, _, _ = self._roundtrip(lines, contig_map={'old_scf1': 'NEW_scf1'}, gene_tag_map={})
        self.assertIn('##sequence-region NEW_scf1 1 1000', content)
        self.assertNotIn('old_scf1', content)

    def test_embedded_fasta_contig_renamed(self):
        lines = ['##gff-version 3\n', '##FASTA\n', '>old_scf1 desc\n', 'ATCG\n']
        content, _, _ = self._roundtrip(lines, contig_map={'old_scf1': 'NEW_scf1'}, gene_tag_map={})
        self.assertIn('>NEW_scf1 desc', content)
        self.assertNotIn('>old_scf1', content)
        self.assertIn('ATCG', content)

    def test_unmapped_values_preserved(self):
        line = 'old_scf1\t.\tCDS\t1\t9\t.\t+\t0\tID=OLD_000001\n'
        content, seqid_changed, attr_renamed = self._roundtrip(
            [line],
            contig_map={'other_scf': 'X'},
            gene_tag_map={'OTHER_000001': 'Y'},
        )
        self.assertIn('old_scf1', content)
        self.assertIn('OLD_000001', content)
        self.assertEqual(seqid_changed, 0)
        self.assertEqual(attr_renamed, 0)

    def test_short_lines_passed_through_unchanged(self):
        """Lines with fewer than 9 tab-separated columns are not modified."""
        line = 'old_scf1\t.\tgene\n'
        content, _, _ = self._roundtrip([line], contig_map={'old_scf1': 'NEW_scf1'}, gene_tag_map={})
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


class TestApplyMapsToGffProkka(TestCase):
    """_apply_maps_to_gff must handle Prokka-specific GFF quirks."""

    def _roundtrip(self, lines, contig_map, gene_tag_map):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as src:
            src.writelines(lines)
            src_path = src.name
        dst_path = src_path + '.out'
        try:
            seqid_changed, attr_renamed = _apply_maps_to_gff(src_path, dst_path, contig_map, gene_tag_map)
            with open(dst_path) as f:
                content = f.read()
            return content, seqid_changed, attr_renamed
        finally:
            os.unlink(src_path)
            if os.path.exists(dst_path):
                os.unlink(dst_path)

    def test_seqid_gnl_prefix_renamed(self):
        """GFF seqid column 'gnl|C|OLD_1' is renamed to the new contig ID."""
        line = 'gnl|C|ALNJDMAK_1\t.\tCDS\t1\t9\t.\t+\t0\tID=OLD_000001;locus_tag=OLD_000001\n'
        content, seqid_changed, attr_renamed = self._roundtrip(
            [line],
            contig_map={'ALNJDMAK_1': 'GENOME_scf1'},
            gene_tag_map={'OLD_000001': 'NEW_000001'},
        )
        self.assertIn('GENOME_scf1\t', content)
        self.assertNotIn('gnl|C|ALNJDMAK_1', content)
        self.assertEqual(seqid_changed, 1)
        self.assertGreater(attr_renamed, 0)

    def test_sequence_region_gnl_prefix_renamed(self):
        """##sequence-region with gnl|C| prefix is renamed."""
        lines = ['##gff-version 3\n', '##sequence-region gnl|C|ALNJDMAK_1 1 40066\n']
        content, _, _ = self._roundtrip(
            lines,
            contig_map={'ALNJDMAK_1': 'GENOME_scf1'},
            gene_tag_map={},
        )
        self.assertIn('##sequence-region GENOME_scf1 1 40066', content)
        self.assertNotIn('gnl|C|ALNJDMAK_1', content)

    def test_embedded_fasta_gnl_prefix_renamed(self):
        """Embedded ##FASTA contig header with gnl|C| prefix is renamed."""
        lines = ['##gff-version 3\n', '##FASTA\n', '>gnl|C|ALNJDMAK_1 desc\n', 'ATCG\n']
        content, _, _ = self._roundtrip(
            lines,
            contig_map={'ALNJDMAK_1': 'GENOME_scf1'},
            gene_tag_map={},
        )
        self.assertIn('>GENOME_scf1 desc', content)
        self.assertNotIn('gnl|C|ALNJDMAK_1', content)

    def test_id_gene_suffix_renamed(self):
        """ID=locus_tag_gene is renamed to ID=new_tag_gene."""
        line = 'scf1\t.\tgene\t1\t9\t.\t+\t.\tID=OLD_000001_gene;locus_tag=OLD_000001\n'
        content, _, _ = self._roundtrip(
            [line],
            contig_map={},
            gene_tag_map={'OLD_000001': 'NEW_000001'},
        )
        self.assertIn('ID=NEW_000001_gene', content)
        self.assertNotIn('ID=OLD_000001_gene', content)

    def test_parent_gene_suffix_renamed(self):
        """Parent=locus_tag_gene is renamed to Parent=new_tag_gene."""
        line = 'scf1\t.\tCDS\t1\t9\t.\t+\t0\tID=OLD_000001;Parent=OLD_000001_gene;locus_tag=OLD_000001\n'
        content, _, _ = self._roundtrip(
            [line],
            contig_map={},
            gene_tag_map={'OLD_000001': 'NEW_000001'},
        )
        self.assertIn('Parent=NEW_000001_gene', content)
        self.assertNotIn('Parent=OLD_000001_gene', content)

    def test_prokka_gene_and_cds_pair_fully_renamed(self):
        """A realistic Prokka gene+CDS pair is fully renamed: seqid, ID, Parent, locus_tag."""
        gene_line = (
            'gnl|C|ALNJDMAK_1\tprokka\tgene\t32\t637\t.\t-\t.'
            '\tID=OLD_00001_gene;locus_tag=OLD_00001\n'
        )
        cds_line = (
            'gnl|C|ALNJDMAK_1\tProdigal:002006\tCDS\t32\t637\t.\t-\t0'
            '\tID=OLD_00001;Parent=OLD_00001_gene;locus_tag=OLD_00001\n'
        )
        contig_map = {'ALNJDMAK_1': 'GENOME_scf1'}
        gene_tag_map = {'OLD_00001': 'NEW_000001'}
        content, _, _ = self._roundtrip([gene_line, cds_line], contig_map, gene_tag_map)

        self.assertNotIn('gnl|C|ALNJDMAK_1', content)
        self.assertIn('GENOME_scf1\t', content)
        self.assertIn('ID=NEW_000001_gene', content)
        self.assertIn('locus_tag=NEW_000001', content)
        self.assertIn('ID=NEW_000001', content)
        self.assertIn('Parent=NEW_000001_gene', content)
        self.assertNotIn('OLD_00001', content)

    def test_pgap_gene_cds_exon_fully_renamed(self):
        """A realistic PGAP gene+CDS+rRNA+exon block is fully renamed."""
        lines = [
            '##gff-version 3\n',
            '##sequence-region NZ_CP015496.1 1 2209387\n',
            'NZ_CP015496.1\tRefSeq\tgene\t1\t301\t.\t+\t.'
            '\tID=gene-RS00005;Name=RS00005;locus_tag=RS00005\n',
            'NZ_CP015496.1\tRefSeq\tCDS\t1\t301\t.\t+\t0'
            '\tID=cds-WP_123.1;Parent=gene-RS00005;locus_tag=RS00005;product=hypothetical\n',
            'NZ_CP015496.1\tRefSeq\tgene\t400\t500\t.\t+\t.'
            '\tID=gene-RS00150;locus_tag=RS00150\n',
            'NZ_CP015496.1\tRefSeq\trRNA\t400\t500\t.\t+\t.'
            '\tID=rna-RS00150;Parent=gene-RS00150;locus_tag=RS00150\n',
            'NZ_CP015496.1\tRefSeq\texon\t400\t500\t.\t+\t.'
            '\tID=exon-RS00150-1;Parent=rna-RS00150\n',
            'NZ_CP015496.1\tRefSeq\tregion\t1\t2209387\t.\t+\t.'
            '\tID=id-NZ_CP015496.1:1..2209387\n',
        ]
        contig_map = {'NZ_CP015496.1': 'GENOME_scf1'}
        gene_tag_map = {'RS00005': 'GENOME_000001', 'RS00150': 'GENOME_000042'}
        content, _, _ = self._roundtrip(lines, contig_map, gene_tag_map)

        # seqid and ##sequence-region renamed
        self.assertNotIn('NZ_CP015496.1\t', content)
        self.assertIn('GENOME_scf1\t', content)
        self.assertIn('##sequence-region GENOME_scf1', content)
        # gene feature IDs renamed
        self.assertIn('ID=gene-GENOME_000001', content)
        self.assertIn('ID=gene-GENOME_000042', content)
        # CDS: protein-accession ID left alone, Parent renamed
        self.assertIn('ID=cds-WP_123.1', content)
        self.assertIn('Parent=gene-GENOME_000001', content)
        # rRNA + exon
        self.assertIn('ID=rna-GENOME_000042', content)
        self.assertIn('Parent=rna-GENOME_000042', content)
        self.assertIn('ID=exon-GENOME_000042-1', content)
        self.assertIn('Parent=gene-GENOME_000042', content)
        # region ID unchanged
        self.assertIn('ID=id-NZ_CP015496.1:1..2209387', content)
        # locus_tag attributes renamed
        self.assertIn('locus_tag=GENOME_000001', content)
        self.assertIn('locus_tag=GENOME_000042', content)
        # no old locus tags remain
        self.assertNotIn('RS00005', content)
        self.assertNotIn('RS00150', content)


# ── _apply_gene_tag_map_to_fasta ──────────────────────────────────────────────

class TestApplyGeneTagMapToFasta(TestCase):
    def _roundtrip(self, content, gene_tag_map):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.faa', delete=False) as src:
            src.write(content)
            src_path = src.name
        dst_path = src_path + '.out'
        try:
            renamed, total = _apply_gene_tag_map_to_fasta(src_path, dst_path, gene_tag_map)
            with open(dst_path) as f:
                result = f.read()
            return result, renamed, total
        finally:
            os.unlink(src_path)
            if os.path.exists(dst_path):
                os.unlink(dst_path)

    def test_basic_rename(self):
        result, renamed, total = self._roundtrip(
            '>OLD_000001 some product\nMPKL\n',
            {'OLD_000001': 'NEW_000001'},
        )
        self.assertEqual(renamed, 1)
        self.assertEqual(total, 1)
        self.assertIn('>NEW_000001 some product', result)
        self.assertNotIn('>OLD_', result)

    def test_description_preserved(self):
        result, _, _ = self._roundtrip(
            '>OLD_000001 hypothetical protein [Organism X]\nMPKL\n',
            {'OLD_000001': 'NEW_000001'},
        )
        self.assertIn('hypothetical protein [Organism X]', result)

    def test_unmapped_tag_unchanged(self):
        result, renamed, total = self._roundtrip('>OLD_000001\nMPKL\n', {'OTHER': 'NEW'})
        self.assertEqual(renamed, 0)
        self.assertEqual(total, 1)
        self.assertIn('>OLD_000001', result)

    def test_sequence_lines_not_modified(self):
        """Sequence lines that look like a locus tag are never touched."""
        result, _, _ = self._roundtrip(
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


# ── _extend_gene_tag_map ──────────────────────────────────────────────────────

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


# ── _resolve_contig_id ────────────────────────────────────────────────────────

class TestResolveContigId(TestCase):
    def test_direct_lookup(self):
        contig_map = {'ALNJDMAK_1': 'NEW_scf1'}
        self.assertEqual(_resolve_contig_id('ALNJDMAK_1', contig_map), 'NEW_scf1')

    def test_gnl_prefix_stripped(self):
        """gnl|C|ALNJDMAK_1 resolves via its bare base ALNJDMAK_1."""
        contig_map = {'ALNJDMAK_1': 'NEW_scf1'}
        self.assertEqual(_resolve_contig_id('gnl|C|ALNJDMAK_1', contig_map), 'NEW_scf1')

    def test_arbitrary_pipe_prefix_stripped(self):
        """Any gnl|X| wrapper is stripped, not just gnl|C|."""
        contig_map = {'ALNJDMAK_1': 'NEW_scf1'}
        self.assertEqual(_resolve_contig_id('gnl|extdb|ALNJDMAK_1', contig_map), 'NEW_scf1')

    def test_returns_none_when_not_found(self):
        self.assertIsNone(_resolve_contig_id('UNKNOWN', {'OTHER': 'X'}))

    def test_returns_none_when_base_not_found(self):
        self.assertIsNone(_resolve_contig_id('gnl|C|UNKNOWN', {'OTHER': 'X'}))


# ── _lookup_gene_tag ──────────────────────────────────────────────────────────

class TestLookupGeneTag(TestCase):
    def test_direct_lookup(self):
        gene_map = {'OLD_000001': 'NEW_000001'}
        self.assertEqual(_lookup_gene_tag('OLD_000001', gene_map), 'NEW_000001')

    def test_gene_suffix(self):
        gene_map = {'OLD_000001': 'NEW_000001'}
        self.assertEqual(_lookup_gene_tag('OLD_000001_gene', gene_map), 'NEW_000001_gene')

    def test_trna_suffix(self):
        gene_map = {'OLD_000042': 'NEW_000042'}
        self.assertEqual(_lookup_gene_tag('OLD_000042_tRNA', gene_map), 'NEW_000042_tRNA')

    def test_rrna_suffix(self):
        gene_map = {'OLD_000007': 'NEW_000007'}
        self.assertEqual(_lookup_gene_tag('OLD_000007_rRNA', gene_map), 'NEW_000007_rRNA')

    def test_repeat_region_suffix(self):
        gene_map = {'OLD_000099': 'NEW_000099'}
        self.assertEqual(_lookup_gene_tag('OLD_000099_repeat_region', gene_map), 'NEW_000099_repeat_region')

    def test_unknown_suffix_returned_unchanged(self):
        gene_map = {'OLD_000001': 'NEW_000001'}
        self.assertEqual(_lookup_gene_tag('OLD_000001_unknownsuffix', gene_map), 'OLD_000001_unknownsuffix')

    def test_not_in_map_returned_unchanged(self):
        self.assertEqual(_lookup_gene_tag('COMPLETELY_OTHER', {}), 'COMPLETELY_OTHER')

    def test_pgap_gene_prefix(self):
        """ID=gene-locus_tag → renamed."""
        gene_map = {'RS00005': 'GENOME_000001'}
        self.assertEqual(_lookup_gene_tag('gene-RS00005', gene_map), 'gene-GENOME_000001')

    def test_pgap_cds_prefix_locus_tag(self):
        """ID=cds-locus_tag (pseudogene) → renamed."""
        gene_map = {'RS00010': 'GENOME_000002'}
        self.assertEqual(_lookup_gene_tag('cds-RS00010', gene_map), 'cds-GENOME_000002')

    def test_pgap_cds_prefix_accession_unchanged(self):
        """ID=cds-WP_accession → left unchanged (not a locus tag)."""
        gene_map = {'RS00005': 'GENOME_000001'}
        self.assertEqual(_lookup_gene_tag('cds-WP_012211193.1', gene_map), 'cds-WP_012211193.1')

    def test_pgap_rna_prefix(self):
        """ID=rna-locus_tag → renamed."""
        gene_map = {'RS00150': 'GENOME_000042'}
        self.assertEqual(_lookup_gene_tag('rna-RS00150', gene_map), 'rna-GENOME_000042')

    def test_pgap_exon_prefix_with_number(self):
        """ID=exon-locus_tag-1 → renamed, exon number preserved."""
        gene_map = {'RS00150': 'GENOME_000042'}
        self.assertEqual(_lookup_gene_tag('exon-RS00150-1', gene_map), 'exon-GENOME_000042-1')

    def test_pgap_exon_higher_number(self):
        """Exon number > 1 is preserved correctly."""
        gene_map = {'RS00150': 'GENOME_000042'}
        self.assertEqual(_lookup_gene_tag('exon-RS00150-3', gene_map), 'exon-GENOME_000042-3')

    def test_pgap_id_region_unchanged(self):
        """ID=id-contig:start..end (region ID, not a locus tag) → left unchanged."""
        gene_map = {'RS00005': 'GENOME_000001'}
        self.assertEqual(_lookup_gene_tag('id-NZ_CP015496.1:1..68603', gene_map),
                         'id-NZ_CP015496.1:1..68603')

    def test_pgap_parent_gene_renamed(self):
        """Parent=gene-locus_tag (CDS pointing to its gene parent) → renamed."""
        gene_map = {'RS00005': 'GENOME_000001'}
        self.assertEqual(_lookup_gene_tag('gene-RS00005', gene_map), 'gene-GENOME_000001')


# ── arx-assigned 5-digit locus tags in annotation files ───────────────────────

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
            renamed, total = _apply_gene_tag_map_to_fasta(src_path, dst_path, lt_map)
            self.assertEqual(renamed, 1)
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


# ── import_genome: Prokka FNA gnl|C| prefix normalisation ─────────────────────

class TestImportGenomeProkkaFnaGnlNorm(TestCase):
    """
    import_genome.py strips the gnl|X| prefix from FNA contig IDs before comparing
    them with GBK contig IDs.  Prokka writes FNA headers as "gnl|C|LOCUS_NAME" but
    the GBK LOCUS line is just "LOCUS_NAME".  This class tests the exact
    normalisation one-liner and the FastaFile.rename_contig_ids round-trip so that
    the fix in import_genome.py remains covered even when real Prokka files are
    unavailable.
    """

    def _write_fna(self, headers: list[str]) -> str:
        """Write a minimal FASTA file and return its path."""
        f = tempfile.NamedTemporaryFile(mode='w', suffix='.fna', delete=False)
        for h in headers:
            f.write(f'>{h} some description\nATCGATCGATCG\n')
        f.close()
        return f.name

    @staticmethod
    def _normalize_ids(fna_ids: list[str]) -> list[str]:
        """The same one-liner used in import_genome.py."""
        return [cid.rsplit('|', 1)[1] if '|' in cid else cid for cid in fna_ids]

    def test_gnl_prefix_stripped_matches_gbk_ids(self):
        """gnl|C| prefixed FNA IDs normalise to the same bare IDs as the GBK."""
        fna_ids = ['gnl|C|ALNJDMAK_1', 'gnl|C|ALNJDMAK_2', 'gnl|C|ALNJDMAK_3']
        gbk_ids = ['ALNJDMAK_1', 'ALNJDMAK_2', 'ALNJDMAK_3']
        self.assertEqual(self._normalize_ids(fna_ids), gbk_ids)

    def test_bare_ids_unchanged(self):
        """Bare IDs (no pipe) are returned as-is — non-Prokka genomes unaffected."""
        ids = ['ALNJDMAK_1', 'ALNJDMAK_2']
        self.assertEqual(self._normalize_ids(ids), ids)

    def test_mixed_ids_normalised(self):
        """Mixed gnl|/bare IDs are all resolved correctly."""
        ids = ['gnl|C|ALNJDMAK_1', 'ALNJDMAK_2']
        self.assertEqual(self._normalize_ids(ids), ['ALNJDMAK_1', 'ALNJDMAK_2'])

    def test_rename_contig_ids_rewrites_gnl_headers(self):
        """
        FastaFile.rename_contig_ids called with the canonical IDs (looked up via
        normalised FNA IDs) produces a FASTA file with plain v3 contig headers.
        """
        fna_path = self._write_fna(['gnl|C|ALNJDMAK_1', 'gnl|C|ALNJDMAK_2'])
        try:
            fna = FastaFile(fna_path)
            fna_ids = fna.get_contig_ids()  # ['gnl|C|ALNJDMAK_1', 'gnl|C|ALNJDMAK_2']
            fna_ids_norm = self._normalize_ids(fna_ids)

            # Simulate what import_genome builds from the GBK (bare IDs in GBK order)
            gbk_ids = ['ALNJDMAK_1', 'ALNJDMAK_2']
            contig_id_map = {gbk_id: f'MYGEN_{i + 1:06d}' for i, gbk_id in enumerate(gbk_ids)}

            out_path = fna_path + '.out'
            fna.rename_contig_ids(
                out=out_path,
                new_ids=[contig_id_map[n] for n in fna_ids_norm],
                update_path=False,
            )

            out_fna = FastaFile(out_path)
            renamed_ids = out_fna.get_contig_ids()
            self.assertEqual(renamed_ids, ['MYGEN_000001', 'MYGEN_000002'])
        finally:
            os.unlink(fna_path)
            if os.path.exists(fna_path + '.out'):
                os.unlink(fna_path + '.out')
