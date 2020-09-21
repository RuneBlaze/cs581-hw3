from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO

import sys
aln_file = sys.argv[1]
out_file = sys.argv[2]

aln = AlignIO.read('Tests/TreeConstruction/msa.phy', 'phylip')