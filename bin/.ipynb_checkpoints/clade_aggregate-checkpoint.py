from pathlib import Path
import pandas as pd
import argparse
import glob
import json


parser = argparse.ArgumentParser(description='Process some files.')
"""
parser.add_argument('files', metavar='file', type=str, nargs='*', default=['*'],
                    help='input file(s), use * to select all files in current directory')

args = parser.parse_args()

if '*' in args.files:
    # Получаем все файлы в текущем каталоге
    files = glob.glob('*')
else:
    files = []
    for pattern in args.files:
        files.extend(glob.glob(pattern))
"""
parser.add_argument("-i", "--input", dest="INPUT", help="The path to the input .json", required=True)
args = parser.parse_args()

INPUT = Path(args.INPUT)
# Теперь files содержит список всех файлов для обработки

#files = [Path(i) for i in files if ".json" in Path(i).name]
files = [i for i in INPUT.glob("*.json")]
print(files)
clades_df = pd.DataFrame()
sample_name_list = []
clades_list = []
for file in files:
    with open(file) as f:
        clade = json.load(f)
    print(file.stem)
    sample_name_list.append(file.stem)
    clades_list.append(clade["final_clade_result"])
clades_df["sample"] = sample_name_list
clades_df["clade"] = clades_list
clades_df.sort_values(by="sample").to_csv("clades_mqc.tsv", sep="\t", index=False)