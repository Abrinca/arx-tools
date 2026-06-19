import os
import re
import warnings

from .utils import GenomeFile, create_replace_function, split_locus_tag

_ID_PREFIX_RE = re.compile(r'([A-Za-z0-9][A-Za-z0-9._-]{3,}_)(\d{5,})')
_STRUCTURAL_PREFIXES = ('cds-', 'gene-', 'rna-', 'exon-', 'id-')
_FREE_TEXT_KEYS = frozenset({'product', 'Note', 'note', 'inference', 'description'})
_MIN_PREFIX_CHARS = 4  # min chars before trailing '_'; filters NCBI 2-letter accession prefixes like WP_, NP_


class NoLocusTagInGffLine(KeyError):
    pass


class GffFile(GenomeFile):
    def rename(
            self,
            out: str,
            new_locus_tag_prefix: str,
            old_locus_tag_prefix: str = None,
            validate: bool = False,
            update_path: bool = True
    ) -> None:
        old_locus_tag_prefix = self._pre_rename_check(out, new_locus_tag_prefix, old_locus_tag_prefix)

        with open(self.path) as in_f:
            content = in_f.read()

        replace_fn = create_replace_function({
            string.format(prefix=old_locus_tag_prefix): string.format(prefix=new_locus_tag_prefix)
            for string in ['-{prefix}', '={prefix}', ':{prefix}', '|{prefix}']
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

    def rename_by_map(self, out: str, lt_map: dict, contig_id_map: dict = None, update_path: bool = True) -> None:
        """Rename locus tags and optionally contig IDs in a GFF file.

        contig_id_map: maps bare contig IDs (without gnl|X| prefix) to new IDs.
        Prokka GFFs use "gnl|X|bare_id" in column 0 and ##sequence-region; the gnl|X|
        prefix is stripped before the lookup so GBK-derived contig maps work directly.
        """
        with open(self.path) as f_in, open(out, 'w') as f_out:
            for line in f_in:
                if line == '##FASTA\n':
                    f_out.write(line)
                    if not contig_id_map:
                        f_out.write(f_in.read())
                        break
                    for line in f_in:
                        if contig_id_map and line.startswith('>'):
                            parts = line[1:].split(None, 1)
                            bare = parts[0].rsplit('|', 1)[1] if '|' in parts[0] else parts[0]
                            if bare in contig_id_map:
                                rest = (' ' + parts[1]) if len(parts) > 1 else '\n'
                                line = f'>{contig_id_map[bare]}{rest}'
                        f_out.write(line)
                    break
                # Rename contig ID in ##sequence-region headers
                if contig_id_map and line.startswith('##sequence-region '):
                    parts = line.split(None, 2)
                    if len(parts) >= 2:
                        old_id = parts[1]
                        bare = old_id.rsplit('|', 1)[1] if '|' in old_id else old_id
                        if bare in contig_id_map:
                            f_out.write(line.replace(old_id, contig_id_map[bare], 1))
                            continue
                if line.startswith('#') or not line.strip():
                    f_out.write(line)
                    continue
                try:
                    data = self._extract_gff_data(line)
                except AssertionError:
                    f_out.write(line)
                    continue
                cols = line.split('\t')
                # Rename contig ID in column 0
                if contig_id_map:
                    old_contig = cols[0]
                    bare = old_contig.rsplit('|', 1)[1] if '|' in old_contig else old_contig
                    if bare in contig_id_map:
                        cols[0] = contig_id_map[bare]
                # Rename locus tag in column 8
                if 'locus_tag' in data:
                    old_tag = data['locus_tag']
                    if old_tag not in lt_map:
                        raise ValueError(f'Locus tag {old_tag!r} not found in lt_map. {self.path=}')
                    cols[8] = cols[8].replace(old_tag, lt_map[old_tag])
                f_out.write('\t'.join(cols))
        if update_path:
            self.path = out

    def get_seqids(self) -> list[str]:
        seen = []
        seen_set = set()
        with open(self.path) as f:
            for line in f:
                if line == '##FASTA\n':
                    break
                if line.startswith('#') or not line.strip():
                    continue
                seqid = line.split('\t', 1)[0]
                if seqid not in seen_set:
                    seen.append(seqid)
                    seen_set.add(seqid)
        return seen

    def find_unexpected_id_prefixes(self, expected_prefix: str) -> dict[str, tuple[int, str]]:
        """Scan every key=value attribute for PREFIX_NNNNN patterns where PREFIX != expected_prefix.

        Skips known free-text keys (product, Note, inference, description) to avoid false positives
        from product descriptions that reference other organisms' gene names. All other keys are
        checked, including old_locus_tag, Alias, and any tool-specific keys.
        Prefixes with fewer than _MIN_PREFIX_CHARS chars before the trailing '_' are skipped
        to avoid flagging short NCBI accession prefixes (WP_, NP_, …) used as CDS feature IDs
        in PGAP-annotated GFFs.
        Returns {unexpected_prefix: (occurrence_count, example_snippet)}.
        """
        counts: dict[str, tuple[int, str]] = {}
        with open(self.path) as f:
            for line in f:
                if line == '##FASTA\n':
                    break
                if line.startswith('#') or not line.strip():
                    continue
                cols = line.split('\t')
                if len(cols) != 9:
                    continue
                for kv in cols[8].split(';'):
                    if '=' not in kv:
                        continue
                    key, _, value = kv.partition('=')
                    if key in _FREE_TEXT_KEYS:
                        continue
                    for m in _ID_PREFIX_RE.finditer(value):
                        prefix = m.group(1)
                        for sp in _STRUCTURAL_PREFIXES:
                            if prefix.startswith(sp):
                                prefix = prefix[len(sp):]
                                break
                        if len(prefix) - 1 < _MIN_PREFIX_CHARS:  # -1 to exclude trailing '_'
                            if prefix != expected_prefix:
                                warnings.warn(
                                    f'Short prefix {prefix!r} found in {self.path!r}; '
                                    f'expected {expected_prefix!r}. This may be an NCBI accession '
                                    f'prefix (e.g. WP_, NP_) or an unexpected short locus tag prefix.',
                                    stacklevel=2,
                                )
                            continue
                        if prefix != expected_prefix:
                            if prefix not in counts:
                                counts[prefix] = (0, f'{key}={value[:m.end()]}')
                            n, ex = counts[prefix]
                            counts[prefix] = (n + 1, ex)
        return counts

    def detect_locus_tag_prefix(self) -> str:
        with open(self.path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                try:
                    locus_tag_prefix, gene_id = self._extract_gff_locus_tag(line)
                except NoLocusTagInGffLine:
                    continue
                return locus_tag_prefix

        raise KeyError(f'Could not extract locus_tag from {self.path=}')

    def validate_locus_tags(self, locus_tag_prefix: str = None):
        with open(self.path) as f:
            for line in f:
                if line == '##FASTA\n':
                    break  # prokka
                if line.startswith('#'):
                    continue
                try:
                    real_locus_tag_prefix, gene_id = self._extract_gff_locus_tag(line)
                except NoLocusTagInGffLine:
                    continue  # in PGAP gffs, some lines contain no locus_tag
                if real_locus_tag_prefix != locus_tag_prefix:
                    raise ValueError(
                        f'Locus tag prefix in {os.path.basename(self.path)!r} does not match: '
                        f'expected {locus_tag_prefix!r}, found {real_locus_tag_prefix!r}. '
                        f'Use rename mode to normalize (--rename on CLI, '
                        f'"Rename locus tags and contig IDs" in web UI).'
                    )

    @staticmethod
    def _extract_gff_data(line: str) -> dict:
        line = line.rstrip('\n').split('\t')
        assert len(line) == 9, f'gff line is malformed! {len(line)=} {line=}'
        return dict(info.split('=', 1) for info in line[8].split(';') if '=' in info)

    @classmethod
    def _extract_gff_locus_tag(cls, line: str) -> (str, str):
        data = cls._extract_gff_data(line)
        if 'locus_tag' not in data:
            raise NoLocusTagInGffLine(f'gff data contains no locus_tag! {line}')
        locus_tag = data['locus_tag']
        locus_tag_prefix, gene_id = split_locus_tag(locus_tag)
        assert ' ' not in locus_tag_prefix, f'The locus_tag may not contain blanks! {locus_tag=}'
        return locus_tag_prefix, gene_id


def rename_gff(
        file: str, out: str,
        new_locus_tag_prefix: str,
        old_locus_tag_prefix: str = None,
        validate: bool = False
):
    """
    Change the locus tags in a general feature format (gff) file

    :param file: input file
    :param out: output file
    :param new_locus_tag_prefix: desired locus tag
    :param old_locus_tag_prefix: locus tag to replace
    :param validate: if true, perform sanity check
    """
    GffFile(
        file=file
    ).rename(
        out=out,
        new_locus_tag_prefix=new_locus_tag_prefix,
        old_locus_tag_prefix=old_locus_tag_prefix,
        validate=validate
    )


def main():
    import fire

    fire.Fire(rename_gff)


if __name__ == '__main__':
    main()
