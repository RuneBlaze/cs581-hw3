from glob import glob
from Bio import AlignIO
import os

inputs = glob("./1000M1/1000M1/*/rose.aln.true.fasta")
inputs += glob("./1000M4/1000M4/*/rose.aln.true.fasta")

for f in inputs:
    print(f"processing: {f}")
    o = os.path.splitext(f)[0] + '.phy'
    with open(f, 'rU') as input_handle, open(o, 'w') as output_handle:
        alignment = AlignIO.parse(input_handle, "fasta")
        AlignIO.write(alignment, output_handle, "phylip")