import json
import os
import re
import warnings
from dataclasses import dataclass, field

from Bio import SeqIO

from .utils import contig_format_to_regex

_LT_DIGITS = 6
_DEFAULT_CONTIG_FORMAT = '_scf{n}'


@dataclass
class V3CheckResult:
    is_v3: bool = True
    has_pending_v3_files: bool = False
    issues: list = field(default_factory=list)
    pending_files: list = field(default_factory=list)

    def summary(self, genome_id: str) -> str:
        if self.has_pending_v3_files:
            names = ', '.join(os.path.basename(f) for f in self.pending_files)
            return f'{genome_id}: PARTIAL UPGRADE — .v3 files exist but not promoted: {names}'
        if self.is_v3:
            return f'{genome_id}: v3 OK'
        return f'{genome_id}: NOT v3 — {"; ".join(self.issues)}'


def check_genome_v3(
    genome_dir: str,
    genome_id: str,
    deep: bool = False,
    contig_format: str = _DEFAULT_CONTIG_FORMAT,
) -> V3CheckResult:
    """
    Check whether a genome folder is v3-compatible.

    Reads genome.json to find file paths; works regardless of subdirectory layout.

    Shallow (default): checks GBK locus_tags + contig IDs and assembly FNA headers.
    Deep: also checks all custom annotation files listed in genome.json.

    Also reports if .v3 intermediate files exist but have not been promoted yet.
    """
    result = V3CheckResult()
    json_path = os.path.join(genome_dir, 'genome.json')

    if not os.path.exists(json_path):
        result.is_v3 = False
        result.issues.append('genome.json not found')
        return result

    with open(json_path) as f:
        genome_json = json.load(f)

    lt_pattern = re.compile(rf'^{re.escape(genome_id)}_\d{{{_LT_DIGITS},}}$')
    contig_pattern = re.compile(rf'^{re.escape(genome_id)}{contig_format_to_regex(contig_format)}$')

    gbk_filename = genome_json.get('cds_tool_gbk_file')
    asm_filename = genome_json.get('assembly_fasta_file')
    custom_annotations = genome_json.get('custom_annotations', [])

    # Detect pending .v3 files
    for filename in ([gbk_filename] if gbk_filename else []) + ([asm_filename] if asm_filename else []):
        v3 = os.path.join(genome_dir, filename) + '.v3'
        if os.path.exists(v3):
            result.has_pending_v3_files = True
            result.pending_files.append(v3)
    if deep:
        for ca in custom_annotations:
            v3 = os.path.join(genome_dir, ca['file']) + '.v3'
            if os.path.exists(v3):
                result.has_pending_v3_files = True
                result.pending_files.append(v3)

    # Check GBK
    if gbk_filename:
        gbk_path = os.path.join(genome_dir, gbk_filename)
        if not os.path.exists(gbk_path):
            result.is_v3 = False
            result.issues.append(f'GBK not found: {gbk_filename}')
        else:
            _check_gbk(gbk_path, lt_pattern, contig_pattern, result)

    # Check assembly FNA
    if asm_filename:
        asm_path = os.path.join(genome_dir, asm_filename)
        if os.path.exists(asm_path):
            _check_fna_headers(asm_path, contig_pattern, result)

    # Deep: check custom annotation files
    if deep:
        for ca in custom_annotations:
            ca_path = os.path.join(genome_dir, ca['file'])
            if os.path.exists(ca_path):
                is_eggnog = ca['type'].startswith('eggnog')
                _check_annotation(ca_path, lt_pattern, is_eggnog, result)

    return result


def _check_gbk(gbk_path: str, lt_pattern, contig_pattern, result: V3CheckResult):
    bad_lts = []
    bad_contigs = []

    with open(gbk_path) as f, warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=UserWarning, module='Bio')
        for rec in SeqIO.parse(f, 'genbank'):
            if not contig_pattern.match(rec.id) and len(bad_contigs) < 3:
                bad_contigs.append(rec.id)
            for feature in rec.features:
                lt = feature.qualifiers.get('locus_tag', [None])[0]
                if lt and not lt_pattern.match(lt) and len(bad_lts) < 3:
                    bad_lts.append(lt)

    if bad_contigs:
        result.is_v3 = False
        result.issues.append(f'GBK contig IDs not v3 (e.g. {", ".join(bad_contigs)})')
    if bad_lts:
        result.is_v3 = False
        result.issues.append(f'GBK locus_tags not v3 (e.g. {", ".join(bad_lts)})')


def _check_fna_headers(fna_path: str, contig_pattern, result: V3CheckResult):
    bad = []
    with open(fna_path) as f:
        for line in f:
            if line.startswith('>'):
                contig_id = line[1:].split(None, 1)[0]
                if not contig_pattern.match(contig_id):
                    bad.append(contig_id)
                    if len(bad) >= 3:
                        break
    if bad:
        result.is_v3 = False
        result.issues.append(f'FNA headers not v3 (e.g. {", ".join(bad)})')


def _check_annotation(ca_path: str, lt_pattern, is_eggnog: bool, result: V3CheckResult):
    bad = []
    with open(ca_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            raw_tag = line.split('\t', 1)[0]
            if is_eggnog and '|' in raw_tag:
                raw_tag = raw_tag.rsplit('|', 1)[1]
            if not lt_pattern.match(raw_tag):
                bad.append(raw_tag)
                if len(bad) >= 3:
                    break
    if bad:
        result.is_v3 = False
        fname = os.path.basename(ca_path)
        result.issues.append(f'{fname}: locus_tags not v3 (e.g. {", ".join(bad)})')
