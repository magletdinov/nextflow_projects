#!/usr/bin/env python3
from pathlib import Path
from argparse import ArgumentParser
import subprocess
import json

#argparse
parser = ArgumentParser()

#argparse
parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="INPUT", help="The path to the input .json", required=True)
parser.add_argument("-o", "--output", dest="OUTPUT", help="The path to the output .json", required=True)
args = parser.parse_args()

INPUT, OUTPUT = Path(args.INPUT), Path(args.OUTPUT)

with open(INPUT) as json_file:
    data = json.load(json_file)
key = [i for i in data.keys()][0]

with open(OUTPUT, 'w') as file:
    json.dump(data[key][8], file)