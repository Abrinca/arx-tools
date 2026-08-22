"""Shared test helpers and fixtures for update_folder_structure tests."""
import json
import os

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

GENOME_ID = 'FAM1079-i1-1.1'


def _write_gbk(path: str, contigs: list[tuple[str, list[str]]]) -> None:
    """Write a minimal GenBank file. contigs: [(contig_id, [locus_tag, ...]), ...]"""
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
    """Write a minimal FASTA file with given contig IDs (with a description field)."""
    with open(path, 'w') as f:
        for contig_id in contigs:
            f.write(f'>{contig_id} some description\nATCGATCGATCG\n')


def _write_annotation(path: str, rows: list[tuple]) -> None:
    """Write a minimal tab-separated annotation file."""
    with open(path, 'w') as f:
        f.write('# comment\n')
        for row in rows:
            f.write('\t'.join(row) + '\n')


def _setup_v2_genome(genome_dir: str) -> None:
    """
    Populate *genome_dir* with a minimal v2 genome (GBK + FNA + KG annotation)
    and write genome.json.  The directory must already exist.
    """
    gbk_path = os.path.join(genome_dir, f'{GENOME_ID}.gbk')
    fna_path = os.path.join(genome_dir, f'{GENOME_ID}.fna')
    ann_path = os.path.join(genome_dir, f'{GENOME_ID}.KG')

    _write_gbk(gbk_path, [
        ('old_contig_1', ['OLD_00001', 'OLD_00002']),
        ('old_contig_2', ['OLD_00003']),
    ])
    _write_fna(fna_path, ['old_contig_1', 'old_contig_2'])
    _write_annotation(ann_path, [
        ('OLD_00001', 'K00001'),
        ('OLD_00002', 'K00002'),
        ('OLD_00003', 'K00003'),
    ])

    genome_json = {
        'identifier': GENOME_ID,
        'cds_tool_gbk_file': f'{GENOME_ID}.gbk',
        'assembly_fasta_file': f'{GENOME_ID}.fna',
        'custom_annotations': [{'file': f'{GENOME_ID}.KG', 'type': 'KG'}],
    }
    with open(os.path.join(genome_dir, 'genome.json'), 'w') as f:
        json.dump(genome_json, f)
