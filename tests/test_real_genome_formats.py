"""
Integration tests against real genome files present on this machine.

Each test is guarded with skipUnless(os.path.exists(...)) so that the test
suite continues to pass even when the external paths are unavailable.

Four distinct annotation-tool output formats appear in the real data:

  Prokka (Abrinca / 20x84EH010_ID1695.003)
    Seqid:    gnl|C|LOCUS_NAME  (prefix stripped by _resolve_contig_id)

  PGAP GCA / WGS (99.0676_GCA_012488935.1)
    Seqid:    AASQEV010000001.1  (NCBI WGS accession, direct match in contig_map)
    Locus ID: gene-B4S48_000001 / cds-EFF5198327.1 / exon-RS-N  (gene-/cds- prefix)

  PGAP GCF / RefSeq (MOD1-EC5047_GCF_002233245.1_ASM223324v1)
    Seqid:    NZ_NLUC01000099.1  (same PGAP code path, NZ_ namespace)
    Locus ID: gene-AS963_RS25430  (RS-style locus tags, same gene- prefix)

  arx in-house / FAM* (FAM20446-i1-1.1)
    Seqid:    FAM20446-i1-1_scf0001  (old 4-digit scf format → rename only contigs)
    Locus ID: gene-FAM20446-i1-1.1_000001  (already 6-digit → identity lt_map)

  FAM1079-i1-1.1 (arx_container + arx_container_perf)
    check_genome_v3 reports NOT v3 for both copies

  Bakta (thomas-1.1)
    Seqid:    contig_1  (plain contig ID, direct match in contig_map)
    Locus ID: thomas-1.1_00005  (bare locus tag, no prefix/suffix → direct lookup)
"""
import os
import re
import tempfile
import unittest

from arx_tools.update_folder_structure import (
    _apply_contig_map_to_fna,
)
from arx_tools.check_v3 import check_genome_v3
from arx_tools.rename_genbank import GenBankFile


# ── path constants ─────────────────────────────────────────────────────────────

_PROKKA_DIR = (
    '/home/sandro/Documents/Abrinca/cases'
    '/20x84EH010_ID1695.003/genomes/20x84EH010_ID1695.003'
)
_PROKKA_GBK = os.path.join(_PROKKA_DIR, '20x84EH010_ID1695.003.gbk')
_PROKKA_FNA = os.path.join(_PROKKA_DIR, '20x84EH010_ID1695.003.fna')
_PROKKA_GFF = os.path.join(_PROKKA_DIR, '20x84EH010_ID1695.003.gff')
_PROKKA_GENOME_ID = '20x84EH010_ID1695.003'

_PGAP_DIR = (
    '/data/repos/arx/arx_container/folder_structure'
    '/organisms/99.0676/genomes/99.0676_GCA_012488935.1'
)
_PGAP_GBK = os.path.join(_PGAP_DIR, '2_cds/99.0676_GCA_012488935.1.gbk')
_PGAP_FNA = os.path.join(_PGAP_DIR, '1_assembly/99.0676_GCA_012488935.1.fna')
_PGAP_GFF = os.path.join(_PGAP_DIR, '2_cds/99.0676_GCA_012488935.1.gff')
_PGAP_GENOME_ID = '99.0676_GCA_012488935.1'

_FAM1079_DIR_CONTAINER = (
    '/data/repos/arx/arx_container/folder_structure'
    '/organisms/FAM1079/genomes/FAM1079-i1-1.1'
)
_FAM1079_DIR_PERF = (
    '/data/repos/arx/arx_container_perf/folder_structure'
    '/organisms/FAM1079/genomes/FAM1079-i1-1.1'
)
_FAM1079_GENOME_ID = 'FAM1079-i1-1.1'

# PGAP GCF / RefSeq :NZ_ contig accessions, RS-style locus tags
_NZ_BASE = '/data/repos/arx/arx_container/folder_structure/organisms'
_NZ_DIR = f'{_NZ_BASE}/MOD1-EC5047/genomes/MOD1-EC5047_GCF_002233245.1_ASM223324v1'
_NZ_GBK = os.path.join(_NZ_DIR, '2_cds/MOD1-EC5047_GCF_002233245.1_ASM223324v1.gbk')
_NZ_GFF = os.path.join(_NZ_DIR, '2_cds/MOD1-EC5047_GCF_002233245.1_ASM223324v1.gff')
_NZ_GENOME_ID = 'MOD1-EC5047_GCF_002233245.1_ASM223324v1'

# arx in-house / FAM* :old 4-digit scf contig IDs, already-6-digit locus tags (identity map)
_FAM20446_DIR = (
    '/data/repos/arx/arx_container/folder_structure'
    '/organisms/FAM20446/genomes/FAM20446-i1-1.1'
)
_FAM20446_GBK = os.path.join(_FAM20446_DIR, '2_cds/FAM20446-i1-1.1.gbk')
_FAM20446_GFF = os.path.join(_FAM20446_DIR, '2_cds/FAM20446-i1-1.1.gff')
_FAM20446_FNA = os.path.join(_FAM20446_DIR, '1_assembly/FAM20446-i1-1.fna')
_FAM20446_GENOME_ID = 'FAM20446-i1-1.1'

# Bakta :plain contig IDs, bare locus tag feature IDs (no Prokka suffix, no PGAP prefix)
_BAKTA_DIR = '/home/sandro/Documents/Abrinca/cases/bakta'
_BAKTA_GBK = os.path.join(_BAKTA_DIR, 'thomas-1.1.gbff')
_BAKTA_FNA = os.path.join(_BAKTA_DIR, 'thomas-1.1.fna')
_BAKTA_GFF = os.path.join(_BAKTA_DIR, 'thomas-1.1.gff3')
_BAKTA_GENOME_ID = 'thomas-1.1'


# ── Prokka tests ───────────────────────────────────────────────────────────────

@unittest.skipUnless(os.path.exists(_PROKKA_GBK), 'Prokka genome not available')
class TestProkkaGnlContigRename(unittest.TestCase):
    """normalize() + _apply_contig_map_to_fna handle Prokka's gnl|C| prefix."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        # Build contig_map from the original GBK (reads only a handful of records)
        out_gbk = os.path.join(self.tmp, 'out.gbk')
        self.contig_map, self.lt_map = GenBankFile(_PROKKA_GBK).normalize(
            out=out_gbk, genome_id=_PROKKA_GENOME_ID
        )
    def tearDown(self):
        self._tmp.cleanup()

    def test_contig_map_keys_are_locus_name(self):
        """contig_map keys must be plain LOCUS names (no gnl| prefix)."""
        for key in self.contig_map:
            self.assertNotIn('|', key, f'contig_map key should not contain pipe: {key!r}')

    def test_contig_map_first_contig_renamed(self):
        """ALNJDMAK_1 → 20x84EH010_ID1695.003_scf1 (first contig)."""
        self.assertIn('ALNJDMAK_1', self.contig_map)
        self.assertEqual(self.contig_map['ALNJDMAK_1'], f'{_PROKKA_GENOME_ID}_scf1')

    def test_fna_gnl_headers_renamed(self):
        """_apply_contig_map_to_fna rewrites gnl|C| prefixed headers to v3 IDs."""
        out_fna = os.path.join(self.tmp, 'out.fna')
        count = _apply_contig_map_to_fna(_PROKKA_FNA, out_fna, self.contig_map)
        self.assertGreater(count, 0, 'Expected at least one contig to be renamed')
        # Verify the first few headers in the output
        with open(out_fna) as f:
            lines = [f.readline() for _ in range(4)]
        headers = [l.strip() for l in lines if l.startswith('>')]
        self.assertTrue(any(_PROKKA_GENOME_ID in h for h in headers),
                        f'No v3 header found in first few lines: {headers}')
        # No gnl| prefix should remain in those lines
        self.assertFalse(any('gnl|' in l for l in lines),
                         f'gnl| prefix still present: {lines}')


# ── PGAP tests ─────────────────────────────────────────────────────────────────

@unittest.skipUnless(os.path.exists(_PGAP_GBK), 'PGAP genome not available')
class TestPgapContigRename(unittest.TestCase):
    """normalize() maps NCBI-accession contig IDs for a PGAP-annotated genome."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        out_gbk = os.path.join(self.tmp, 'out.gbk')
        self.contig_map, self.lt_map = GenBankFile(_PGAP_GBK).normalize(
            out=out_gbk, genome_id=_PGAP_GENOME_ID
        )

    def tearDown(self):
        self._tmp.cleanup()

    def test_contig_map_first_contig(self):
        """AASQEV010000001.1 → 99.0676_GCA_012488935.1_scf1."""
        self.assertIn('AASQEV010000001.1', self.contig_map)
        self.assertEqual(self.contig_map['AASQEV010000001.1'], f'{_PGAP_GENOME_ID}_scf1')

    def test_lt_map_uses_6_digit_values(self):
        """All values in lt_map must be 6-digit v3 locus tags."""
        import re
        v3_pattern = re.compile(rf'^{re.escape(_PGAP_GENOME_ID)}_\d{{6}}$')
        for new_lt in self.lt_map.values():
            self.assertRegex(new_lt, v3_pattern)


# ── FAM1079 v3 check ───────────────────────────────────────────────────────────

@unittest.skipUnless(os.path.exists(_FAM1079_DIR_CONTAINER), 'FAM1079 arx_container not available')
class TestFam1079NotV3Container(unittest.TestCase):
    """FAM1079-i1-1.1 in arx_container is NOT v3 (old contig IDs, old locus tags)."""

    def test_check_genome_v3_reports_not_v3(self):
        result = check_genome_v3(_FAM1079_DIR_CONTAINER, _FAM1079_GENOME_ID)
        self.assertFalse(result.is_v3,
                         f'Expected NOT v3 but got is_v3=True; issues: {result.issues}')

    def test_check_reports_contig_or_locus_issue(self):
        result = check_genome_v3(_FAM1079_DIR_CONTAINER, _FAM1079_GENOME_ID)
        issue_text = ' '.join(result.issues).lower()
        self.assertTrue(
            'contig' in issue_text or 'locus' in issue_text,
            f'Expected a contig/locus issue but got: {result.issues}',
        )


@unittest.skipUnless(os.path.exists(_FAM1079_DIR_PERF), 'FAM1079 arx_container_perf not available')
class TestFam1079NotV3Perf(unittest.TestCase):
    """FAM1079-i1-1.1 in arx_container_perf is also NOT v3."""

    def test_check_genome_v3_reports_not_v3(self):
        result = check_genome_v3(_FAM1079_DIR_PERF, _FAM1079_GENOME_ID)
        self.assertFalse(result.is_v3,
                         f'Expected NOT v3 but got is_v3=True; issues: {result.issues}')


# ── PGAP GCF / RefSeq (NZ_ accessions, RS-style locus tags) ───────────────────

@unittest.skipUnless(os.path.exists(_NZ_GBK), 'NZ_ RefSeq genome not available')
class TestPgapNZContigRename(unittest.TestCase):
    """
    PGAP GCF genomes use NZ_ accession seqids and RS-style locus tags.
    The GFF ID scheme (gene-/cds-) is identical to GCA; this test confirms the
    same code path works with NZ_ prefixes.
    """

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        out_gbk = os.path.join(self.tmp, 'out.gbk')
        self.contig_map, self.lt_map = GenBankFile(_NZ_GBK).normalize(
            out=out_gbk, genome_id=_NZ_GENOME_ID
        )

    def tearDown(self):
        self._tmp.cleanup()

    def test_contig_map_first_contig_uses_nz_key(self):
        """contig_map keys are NZ_ accessions; values are v3 scf IDs."""
        nz_keys = [k for k in self.contig_map if k.startswith('NZ_')]
        self.assertGreater(len(nz_keys), 0, 'Expected NZ_ accession keys in contig_map')
        first_nz = nz_keys[0]
        self.assertTrue(self.contig_map[first_nz].startswith(_NZ_GENOME_ID),
                        f'v3 contig ID should start with genome_id: {self.contig_map[first_nz]}')

    def test_lt_map_uses_6_digit_values(self):
        """All values in lt_map must be 6-digit v3 locus tags."""
        import re
        v3_pattern = re.compile(rf'^{re.escape(_NZ_GENOME_ID)}_\d{{6}}$')
        for new_lt in self.lt_map.values():
            self.assertRegex(new_lt, v3_pattern)


# ── arx in-house FAM* (identity lt_map, only contig rename) ───────────────────

@unittest.skipUnless(os.path.exists(_FAM20446_GBK), 'FAM20446 genome not available')
class TestFamInHouseContigRename(unittest.TestCase):
    """
    arx in-house FAM* genomes already carry 6-digit locus tags; normalize()
    produces an identity lt_map.  Only contig IDs need renaming (scf0001 → scf1).
    """

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        out_gbk = os.path.join(self.tmp, 'out.gbk')
        self.contig_map, self.lt_map = GenBankFile(_FAM20446_GBK).normalize(
            out=out_gbk, genome_id=_FAM20446_GENOME_ID
        )

    def tearDown(self):
        self._tmp.cleanup()

    def test_lt_map_is_identity(self):
        """Locus tags are already v3; every old tag maps to itself."""
        self.assertTrue(all(k == v for k, v in self.lt_map.items()),
                        'Expected identity lt_map for FAM20446')

    def test_contig_map_renames_scf0001_to_scf1(self):
        """Old 4-digit scf IDs are renamed to 1-digit v3 scf IDs."""
        old_keys = [k for k in self.contig_map if 'scf0' in k]
        self.assertGreater(len(old_keys), 0, 'Expected scf0NNN keys in contig_map')
        for old, new in self.contig_map.items():
            self.assertIn('scf0', old, f'Expected old 4-digit scf key: {old}')
            self.assertRegex(new, rf'^{re.escape(_FAM20446_GENOME_ID)}_scf\d+$')


# ── Bakta (plain contig IDs, bare locus tag feature IDs) ──────────────────────

@unittest.skipUnless(os.path.exists(_BAKTA_GBK), 'Bakta genome not available')
class TestBaktaContigRename(unittest.TestCase):
    """Bakta GBK uses plain contig_N IDs; normalize() maps them to v3 scf IDs."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        out_gbk = os.path.join(self.tmp, 'out.gbk')
        self.contig_map, self.lt_map = GenBankFile(_BAKTA_GBK).normalize(
            out=out_gbk, genome_id=_BAKTA_GENOME_ID
        )

    def tearDown(self):
        self._tmp.cleanup()

    def test_contig_map_first_contig(self):
        """contig_1 → thomas-1.1_scf1."""
        self.assertIn('contig_1', self.contig_map)
        self.assertEqual(self.contig_map['contig_1'], f'{_BAKTA_GENOME_ID}_scf1')

    def test_lt_map_uses_6_digit_values(self):
        """5-digit Bakta locus tags are mapped to 6-digit v3 locus tags."""
        v3_pattern = re.compile(rf'^{re.escape(_BAKTA_GENOME_ID)}_\d{{6}}$')
        for new_lt in self.lt_map.values():
            self.assertRegex(new_lt, v3_pattern)


@unittest.skipUnless(os.path.exists(_BAKTA_FNA), 'Bakta FNA not available')
class TestBaktaFnaRename(unittest.TestCase):
    """_apply_contig_map_to_fna renames plain contig_N headers to v3 scf IDs."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        out_gbk = os.path.join(self.tmp, 'out.gbk')
        self.contig_map, _ = GenBankFile(_BAKTA_GBK).normalize(
            out=out_gbk, genome_id=_BAKTA_GENOME_ID
        )

    def tearDown(self):
        self._tmp.cleanup()

    def test_fna_headers_renamed(self):
        """contig_1 header is replaced with thomas-1.1_scf1."""
        out_fna = os.path.join(self.tmp, 'out.fna')
        count = _apply_contig_map_to_fna(_BAKTA_FNA, out_fna, self.contig_map)
        self.assertGreater(count, 0, 'Expected at least one contig to be renamed')
        with open(out_fna) as f:
            first_header = f.readline().strip()
        self.assertTrue(first_header.startswith(f'>{_BAKTA_GENOME_ID}_scf'),
                        f'Unexpected FNA header: {first_header!r}')
        self.assertNotIn('contig_', first_header)


