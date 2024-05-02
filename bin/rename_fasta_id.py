#!/usr/bin/env python3

from Bio import SeqIO
from pathlib import Path
from argparse import ArgumentParser

#argparse
parser = ArgumentParser()

parser.add_argument("-i", "--input", dest="INPUT", help="The path to the input fasta file", required=True)
parser.add_argument("-o", "--output", dest="OUTPUT", help="The path to the output new fasta file", required=True)

args = parser.parse_args()

INPUT = Path(args.INPUT)
OUTPUT = Path(args.OUTPUT)

#1-base sistem
interesting_coord = {
    "locus_on_ORF21":(33658, 33812),
    "locus_on_ORF22":(37829, 38148),
    "locus_on_ORF29":(52307, 52489),
    "locus_on_ORF38":(69268, 69516),
    "locus_on_ORF55":(98375, 98549),
    "locus_on_ORF67":(114552, 114719)}

fa = []
for seq_record in SeqIO.parse(INPUT, "fasta"):
    start = interesting_coord[seq_record.id][0]
    stop = interesting_coord[seq_record.id][1]
    seq_record.description = f"{start}:{stop}_in_NC_001348.1"
    fa.append(seq_record)
SeqIO.write(fa, OUTPUT, "fasta")