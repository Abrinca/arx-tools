import os
import re
import logging

from Bio import SeqIO, SeqRecord, SeqFeature
from .utils import GenomeFile, query_int, entrez_organism_to_taxid, date_to_string, datetime, create_replace_function, \
    split_locus_tag, contig_format_to_regex
from .genbank_to_fasta import GenBankToFasta

# ── GenBank header text-patching helpers ──────────────────────────────────────
# BioPython round-trips lose LOCUS metadata (topology, division, date) from
# Prokka multi-line LOCUS headers, and drop secondary ACCESSION entries due to
# a key mismatch (parser writes 'accessions', writer reads 'accession').
#
# normalize_contigs uses pure text replacement (no BioPython write at all).
# normalize uses BioPython write for feature edits then text-patches headers.

_DATE_RE = re.compile(
    r'\b(\d{2}-(?:JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)-\d{4})\b'
)
_TOPOLOGY_RE = re.compile(r'\b(circular|linear)\b')
_DIVISION_RE = re.compile(
    r'\b(BCT|VRL|PHG|SYN|UNA|RNA|CON|ROD|MAM|INV|PLN|PRI|HUM|ENV|PAT|'
    r'EST|STS|GSS|HTC|TSA|WGS|PRO|UNK)\b'
)
_KNOWN_TOPOLOGY = frozenset({'circular', 'linear'})
_KNOWN_DIVISION = frozenset({
    'BCT', 'VRL', 'PHG', 'SYN', 'UNA', 'RNA', 'CON', 'ROD',
    'MAM', 'INV', 'PLN', 'PRI', 'HUM', 'ENV', 'PAT', 'EST',
    'STS', 'GSS', 'HTC', 'TSA', 'WGS', 'PRO', 'UNK',
})

# Matches a standard or Prokka LOCUS line: name, length, unit, rest-of-line
_LOCUS_LINE_RE = re.compile(r'^LOCUS\s+(\S+)\s+(\d+)\s+(bp|aa)\s+(.*?)$')


def _patch_contig_id_lines(lines: list[str], name_map: dict[str, str]) -> list[str]:
    """
    Rename contig IDs in GenBank file lines by pure text replacement.

    Replaces the LOCUS name field, the primary ACCESSION entry, and the VERSION
    entry.  Handles Prokka multi-line LOCUS (date continuation on next line).

    name_map maps old rec.name (= LOCUS name token) to the new contig ID.
    """
    out: list[str] = []
    i = 0
    current_new_id: str | None = None

    while i < len(lines):
        line = lines[i]

        if line.startswith('LOCUS'):
            m = _LOCUS_LINE_RE.match(line.rstrip('\n'))
            if m:
                old_name, length, unit, rest = m.groups()
                current_new_id = name_map.get(old_name, old_name)

                # Absorb a Prokka date-continuation line (date wrapped to next line)
                if (i + 1 < len(lines)
                        and lines[i + 1].startswith('            ')
                        and _DATE_RE.search(lines[i + 1])):
                    rest = rest.rstrip() + ' ' + lines[i + 1].strip()
                    i += 1

                # BioPython uses a 28-char name_length field: name + right-justified length.
                # If the new name is long, extend the field so length always has a space before it.
                field_width = max(28, len(current_new_id) + 1 + len(length))
                name_length = length.rjust(field_width)
                name_length = current_new_id + name_length[len(current_new_id):]

                out.append(f'LOCUS       {name_length} {unit}    {rest.rstrip()}\n')
            else:
                out.append(line)

        elif line.startswith('ACCESSION'):
            acc_rest = line[len('ACCESSION'):].strip()
            tokens = acc_rest.split()
            if tokens:
                new_primary = name_map.get(tokens[0], current_new_id or tokens[0])
                out.append('ACCESSION   ' + ' '.join([new_primary] + tokens[1:]) + '\n')
            elif current_new_id:
                out.append(f'ACCESSION   {current_new_id}\n')
            else:
                out.append(line)

        elif line.startswith('VERSION'):
            out.append(f'VERSION     {current_new_id}\n' if current_new_id else line)

        elif line.startswith('//'):
            out.append(line)
            current_new_id = None

        else:
            out.append(line)

        i += 1

    return out


def _read_accession_lines(path: str) -> list[str]:
    """Return the raw ACCESSION line (no trailing newline) for each record."""
    results: list[str] = []
    current: str | None = None

    with open(path) as fh:
        for line in fh:
            if line.startswith('LOCUS'):
                current = None
            elif line.startswith('ACCESSION') and current is None:
                current = line.rstrip('\n')
            elif line.startswith('//'):
                results.append(current or '')
                current = None

    return results


def _patch_accession_lines(
    lines: list[str],
    original_accession_lines: list[str],
    name_map: dict[str, str],
) -> list[str]:
    """
    In BioPython-written GBK lines, replace each ACCESSION primary with the new
    contig ID and restore any secondary accessions from the original file.

    name_map maps old rec.name (original primary accession) to the new contig ID.
    """
    out: list[str] = []
    rec_index = -1

    for line in lines:
        if line.startswith('LOCUS'):
            rec_index += 1
            out.append(line)
        elif line.startswith('ACCESSION') and 0 <= rec_index < len(original_accession_lines):
            orig = original_accession_lines[rec_index]
            if orig:
                orig_tokens = orig[len('ACCESSION'):].strip().split()
                if orig_tokens:
                    new_primary = name_map.get(orig_tokens[0], '')
                    if not new_primary:
                        bio_tokens = line[len('ACCESSION'):].strip().split()
                        new_primary = bio_tokens[0] if bio_tokens else ''
                    out.append('ACCESSION   ' + ' '.join([new_primary] + orig_tokens[1:]) + '\n')
                else:
                    out.append(line)
            else:
                out.append(line)
        else:
            out.append(line)

    return out


def _parse_locus_fields_raw(gbk_path: str) -> list[dict]:
    """
    Text-scan a GenBank file and return one dict per record with the raw LOCUS
    line fields (topology, data_file_division, date), handling the Prokka
    multi-line LOCUS format where the date wraps onto the next line.
    Used by normalize() to restore LOCUS metadata after BioPython write.
    """
    results = []
    current: dict | None = None
    in_locus_continuation = False

    with open(gbk_path) as fh:
        for line in fh:
            if line.startswith('LOCUS'):
                current = {'topology': None, 'data_file_division': None, 'date': None}
                in_locus_continuation = True
                _scan_locus_tokens(line, current)
            elif in_locus_continuation and line[:1] == ' ':
                _scan_locus_tokens(line, current)
                in_locus_continuation = False
            else:
                in_locus_continuation = False
                if line.startswith('//') and current is not None:
                    results.append(current)
                    current = None

    if current is not None:
        results.append(current)
    return results


def _scan_locus_tokens(line: str, acc: dict) -> None:
    """Fill missing topology/division/date slots in acc from tokens in line."""
    if acc['topology'] is None:
        m = _TOPOLOGY_RE.search(line)
        if m:
            acc['topology'] = m.group(1)
    if acc['data_file_division'] is None:
        m = _DIVISION_RE.search(line)
        if m:
            acc['data_file_division'] = m.group(1)
    if acc['date'] is None:
        m = _DATE_RE.search(line)
        if m:
            acc['date'] = m.group(1)


def _restore_locus_fields(rec, raw: dict) -> None:
    """
    Apply raw-extracted LOCUS fields to a BioPython SeqRecord, correcting
    misparses from multi-line LOCUS headers.  Used by normalize().
    """
    if raw['topology'] and not rec.annotations.get('topology'):
        rec.annotations['topology'] = raw['topology']
    if rec.annotations.get('data_file_division') in _KNOWN_TOPOLOGY and raw['data_file_division']:
        rec.annotations['data_file_division'] = raw['data_file_division']
    date = rec.annotations.get('date', '')
    if (not date or date in _KNOWN_DIVISION or date == '01-JAN-1980') and raw['date']:
        rec.annotations['date'] = raw['date']


class GenBankFile(GenomeFile):
    def rename(
            self,
            out: str,
            new_locus_tag_prefix: str,
            old_locus_tag_prefix: str = None,
            validate: bool = False,
            scf_prefix: str = None,
            scf_leading_zeroes: int = None,
            update_path: bool = True
    ) -> None:
        old_locus_tag_prefix = self._pre_rename_check(out, new_locus_tag_prefix, old_locus_tag_prefix)

        with open(self.path) as in_f:
            content = in_f.readlines()

        if scf_prefix:
            # change scf prefixes
            if type(scf_leading_zeroes) is int and scf_leading_zeroes > 1:
                format = lambda c: f'VERSION     {scf_prefix}{str(c).zfill(scf_leading_zeroes)}\n'
            else:
                format = lambda c: f'VERSION     {scf_prefix}{c}\n'

            counter = 0
            for i, line in enumerate(content):
                if line.startswith('VERSION'):
                    counter += 1
                    content[i] = format(counter)

        content = ''.join(content)
        # change locus tags
        replace_fn = create_replace_function({
            string.format(prefix=old_locus_tag_prefix): string.format(prefix=new_locus_tag_prefix)
            for string in ['/locus_tag="{prefix}', '/protein_id="extdb:{prefix}', ':{prefix}']
        })

        old_hash = hash(content)
        content = replace_fn(content)
        assert old_hash != hash(content), f'The content of {self.path=} has not changed!'

        assert new_locus_tag_prefix in content, f'Something went wrong: did not replace anything!'

        with open(out, 'w') as out_f:
            out_f.write(content)

        if update_path:
            self.path = out

        if validate:
            self.validate_locus_tags(locus_tag_prefix=new_locus_tag_prefix)

    def create_fna(self, fna: str):
        GenBankToFasta.create_fna(gbk=self.path, out=fna)

    def create_gff(self, gff: str):
        GenBankToFasta.create_gff(gbk=self.path, out=gff)

    def create_ffn(self, ffn: str):
        GenBankToFasta.convert(gbk=self.path, out=ffn, format='ffn')

    def create_faa(self, faa: str):
        GenBankToFasta.convert(gbk=self.path, out=faa, format='faa')

    def normalize_contigs(self, out: str, genome_id: str,
                              contig_ids: list[str] = None,
                              contig_format: str = '_scf{n}') -> dict:
        """
        Canonicalize only contig IDs in a GenBank file using pure text replacement.
        Locus tags are left untouched. All LOCUS metadata is preserved exactly.

        Returns contig_map {old_id: new_id}.
        """
        with open(self.path) as f:
            all_records = list(SeqIO.parse(f, 'genbank'))

        if contig_ids is not None:
            if len(contig_ids) != len(all_records):
                raise ValueError(
                    f'contig_ids count ({len(contig_ids)}) does not match number of records '
                    f'({len(all_records)}) in {self.path}')

        name_map: dict[str, str] = {}
        contig_map: dict[str, str] = {}
        for i, rec in enumerate(all_records):
            new_id = contig_ids[i] if contig_ids is not None else f'{genome_id}{contig_format.format(n=i + 1)}'
            name_map[rec.name] = new_id
            contig_map[rec.id] = new_id

        with open(self.path) as f:
            lines = f.readlines()

        with open(out, 'w') as f:
            f.writelines(_patch_contig_id_lines(lines, name_map))

        return contig_map

    def normalize(self, out: str, genome_id: str,
                  contig_ids: list[str] = None,
                  contig_format: str = '_scf{n}') -> tuple[dict, dict]:
        """
        Canonicalize contig IDs and locus_tags in a GenBank file.

        Contigs are renamed using contig_ids (if provided, matched by position) or
        auto-generated as {genome_id}{contig_format.format(n=counter)}.
        Locus tags are renamed to {genome_id}_{n:06d}.
        Old identifiers are preserved as old_locus_tag qualifiers.

        Returns (contig_map, lt_map) where both are {old_id: new_id} dicts.
        """
        raw_locus = _parse_locus_fields_raw(self.path)
        original_accession_lines = _read_accession_lines(self.path)
        lt_counter = 0
        lt_map = {}
        contig_map = {}
        name_map: dict[str, str] = {}
        records = []

        with open(self.path) as f:
            all_records = list(SeqIO.parse(f, 'genbank'))

        if contig_ids is not None:
            if len(contig_ids) != len(all_records):
                raise ValueError(
                    f'contig_ids count ({len(contig_ids)}) does not match number of records '
                    f'({len(all_records)}) in {self.path}')

        for i, rec in enumerate(all_records):
            if contig_ids is not None:
                new_contig_id = contig_ids[i]
            else:
                new_contig_id = f'{genome_id}{contig_format.format(n=i + 1)}'
            contig_map[rec.id] = new_contig_id
            name_map[rec.name] = new_contig_id
            rec.id = new_contig_id
            rec.name = new_contig_id[:16]  # LOCUS line is limited to 16 chars in genbank format
            if i < len(raw_locus):
                _restore_locus_fields(rec, raw_locus[i])

            for feature in rec.features:
                if 'locus_tag' not in feature.qualifiers:
                    continue
                old_lt = feature.qualifiers['locus_tag'][0]
                if old_lt not in lt_map:
                    lt_counter += 1
                    lt_map[old_lt] = f'{genome_id}_{str(lt_counter).zfill(6)}'
                new_lt = lt_map[old_lt]
                feature.qualifiers['locus_tag'] = [new_lt]
                if 'old_locus_tag' not in feature.qualifiers:
                    feature.qualifiers['old_locus_tag'] = [old_lt]
                elif old_lt not in feature.qualifiers['old_locus_tag']:
                    feature.qualifiers['old_locus_tag'].append(old_lt)
                if feature.qualifiers.get('protein_id', [None])[0] == f'C:{old_lt}':
                    feature.qualifiers['protein_id'] = [f'C:{new_lt}']

            records.append(rec)

        import io
        buf = io.StringIO()
        SeqIO.write(records, buf, 'genbank')
        bio_lines = buf.getvalue().splitlines(keepends=True)
        patched = _patch_accession_lines(bio_lines, original_accession_lines, name_map)

        with open(out, 'w') as f:
            f.writelines(patched)

        return contig_map, lt_map

    def get_contig_ids(self) -> list[str]:
        with open(self.path) as f:
            return [rec.id for rec in SeqIO.parse(f, 'genbank')]

    def validate_contig_ids(self, genome_id: str, contig_format: str = '_scf{n}') -> None:
        """Check that all contig IDs match {genome_id}{contig_format}. Raise ValueError if not."""
        pattern = re.compile(rf'^{re.escape(genome_id)}{contig_format_to_regex(contig_format)}$')
        with open(self.path) as f:
            for rec in SeqIO.parse(f, 'genbank'):
                if not pattern.match(rec.id):
                    raise ValueError(
                        f'Contig ID {rec.id!r} in {os.path.basename(self.path)!r} does not match '
                        f'expected format {genome_id!r} + {contig_format!r} '
                        f'(e.g. {genome_id!r}_scf1). '
                        f'Use rename mode to normalize (--rename on CLI, '
                        f'"Rename locus tags and contig IDs" in web UI).'
                    )

    def validate_locus_tags(self, locus_tag_prefix: str = None):
        if locus_tag_prefix is None:
            locus_tag_prefix = self.detect_locus_tag_prefix()

        with open(self.path) as f:
            for rec in SeqIO.parse(f, "genbank"):
                for feature in rec.features:
                    locus_tag = feature.qualifiers.get('locus_tag')
                    if locus_tag is not None:
                        locus_tag = locus_tag[0]
                        real_locus_tag_prefix, gene_id = split_locus_tag(locus_tag)
                        if real_locus_tag_prefix != locus_tag_prefix:
                            raise ValueError(
                                f'Locus tag prefix in {os.path.basename(self.path)!r} does not match: '
                                f'expected {locus_tag_prefix!r}, found {real_locus_tag_prefix!r}. '
                                f'Use rename mode to normalize (--rename on CLI, '
                                f'"Rename locus tags and contig IDs" in web UI).'
                            )
                        if not gene_id.isdigit():
                            raise ValueError(
                                f'Malformed locus tag {locus_tag!r} in {os.path.basename(self.path)!r}: '
                                f'expected format {locus_tag_prefix!r}_[0-9]+.'
                            )

    def metadata(self) -> (dict, dict):
        organism_data, genome_data = {}, {}

        try:
            organism_data['taxid'] = self.taxid(raise_error=True)
        except Exception as e:
            logging.warning(
                f'Could not determine taxid from {os.path.basename(self.path)}: {e}. '
                f'Provide it manually via organism.json or the taxid form field.'
            )

        rec, feature = self._get_first_gbk_rec_feature(gbk=self.path)
        comment = rec.annotations.get('comment', '')
        if 'prokka' in comment:
            genome_data.update(
                cds_tool='prokka',
                cds_tool_date=date_to_string(datetime.strptime(rec.annotations['date'], '%d-%b-%Y')),
                cds_tool_version=comment.split('prokka ')[-1].split(' ')[0]
            )
        elif 'PGAP' in comment:
            pgap_comment = rec.annotations['structured_comment']['Genome-Annotation-Data']
            genome_data.update(
                cds_tool='PGAP',
                cds_tool_date=date_to_string(datetime.strptime(pgap_comment['Annotation Date'], '%m/%d/%Y %H:%M:%S')),
                cds_tool_version=pgap_comment['Annotation Software revision']
            )
        elif 'Bakta' in comment:
            bakta_comment = rec.annotations['comment']
            bakta_version = bakta_comment.split('Software: ',1)[1].split('\n',1)[0]
            bakta_db_verson = bakta_comment.split('Database: ',1)[1].split('\n',1)[0].replace(', ', '')
            genome_data.update(
                cds_tool='Bakta',
                cds_tool_date=date_to_string(datetime.strptime(rec.annotations['date'], '%d-%b-%Y')),  #08-SEP-2023
                cds_tool_version=f'{bakta_version}:{bakta_db_verson}'
            )
        else:
            logging.warning(f'Failed to discover annotation information from {self.path=}')

        # get  BioProject / BioSample metadata
        if hasattr(rec, 'dbxrefs'):
            for dbxref in rec.dbxrefs:
                if dbxref.startswith('BioProject:'):
                    genome_data['bioproject_accession'] = dbxref.removeprefix('BioProject:')
                if dbxref.startswith('BioSample:'):
                    genome_data['biosample_accession'] = dbxref.removeprefix('BioSample:')

        return organism_data, genome_data

    def taxid(self, raise_error=True, sample_name: str = None) -> int:
        rec, feature = self._get_first_gbk_rec_feature(gbk=self.path)

        # PGAP files contain taxid in first feature (db_xref qualifier)
        db_xref = feature.qualifiers.get('db_xref')
        if type(db_xref) is list and len(db_xref) == 1:
            taxid = db_xref[0].removeprefix('taxon:')
            assert taxid.isdigit(), f'Failed to extract taxid from db_xref in {self.path=}'
            return int(taxid)

        # prokka files contain organism name, which can be turned into taxid using Entrez
        organism = feature.qualifiers.get('organism')
        if type(organism) is list and len(organism) == 1:
            if organism[0] == 'Genus species':
                raise AssertionError(
                    f'GBK organism is the Prokka placeholder "Genus species": Prokka was run without '
                    f'--genus/--species so the taxid cannot be auto-detected. '
                    f'Upload an organism.json alongside the GBK with {{"taxid": <int>}}. '
                    f'Find your taxid at https://www.ncbi.nlm.nih.gov/taxonomy'
                )
            return entrez_organism_to_taxid(organism[0])

        # ask for input or raise error
        if not raise_error:
            taxid = query_int(
                question=f'What is the taxid of {sample_name if sample_name else self.path}?',
                error_msg=f'The response must be an integer. Hint: taxid=1427524 is "mixed sample".'
            )
            assert taxid != 0, f'0 is not a valid taxid.'
            return taxid
        else:
            raise AssertionError(f'Failed to extract taxid from {self.path=}')

    def detect_locus_tag_prefix(self) -> str:
        strain, locus_tag_prefix = self.detect_strain_locus_tag_prefix()
        return locus_tag_prefix

    def detect_strain_locus_tag_prefix(self) -> (str, str):
        strain, locus_tag = None, None
        with open(self.path) as f:
            for rec in SeqIO.parse(f, "genbank"):
                for feature in rec.features:
                    if strain is None:
                        strain = feature.qualifiers.get('strain')
                    if locus_tag is None:
                        locus_tag = feature.qualifiers.get('locus_tag')
                    if strain is not None and locus_tag is not None:
                        break

        if locus_tag is None:
            raise ValueError(
                f'Could not read genome name from {self.path!r}: '
                f'no /locus_tag= qualifier found in any feature.'
            )

        locus_tag_prefix, gene_id = split_locus_tag(locus_tag[0])

        if type(strain) is list and len(strain) == 1 and type(strain[0]) is str:
            strain = strain[0]
        else:
            strain = os.environ.get('STRAIN', None)
            if strain is None:
                strain = locus_tag_prefix.rstrip('_')
                logging.info(
                    f'No /strain= qualifier in {self.path!r}, '
                    f'falling back to locus tag prefix as organism name: {strain!r}'
                )

        return strain, locus_tag_prefix

    @staticmethod
    def _get_first_gbk_rec_feature(gbk: str) -> (SeqRecord, SeqFeature):
        with open(gbk) as f:
            for rec in SeqIO.parse(f, "genbank"):
                for feature in rec.features:
                    return rec, feature
        raise AssertionError(f'Failed to get rec and feature from {gbk=}')


def rename_genbank(
        file: str, out: str,
        new_locus_tag_prefix: str, old_locus_tag_prefix: str = None, validate: bool = False,
        scf_prefix: str = None, scf_leading_zeroes: int = None,
):
    """
    Change the locus tags in a GenBank file

    :param file: input file
    :param out: output file
    :param new_locus_tag_prefix: desired locus tag
    :param old_locus_tag_prefix: locus tag to replace
    :param validate: if true, perform sanity check on locus tags
    :param scf_prefix: desired scaffold prefix (optional)
    :param scf_leading_zeroes: format scaffold counter with leading zeroes. e.g.: 5 -> >PREFIX_00001 (optional)
    """

    GenBankFile(
        file=file
    ).rename(
        out=out,
        new_locus_tag_prefix=new_locus_tag_prefix,
        old_locus_tag_prefix=old_locus_tag_prefix,
        validate=validate,
        scf_prefix=scf_prefix,
        scf_leading_zeroes=scf_leading_zeroes
    )


def main():
    import fire

    fire.Fire(rename_genbank)


if __name__ == '__main__':
    main()
