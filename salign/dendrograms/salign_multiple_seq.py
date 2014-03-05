# Illustrates the SALIGN multiple sequence alignment
import sys

from modeller import *

log.verbose()
env = environ()
env.io.atom_files_directory = ['.', '../../pdbs']

aln = alignment(env, file=sys.argv[1], alignment_format="FASTA")

aln.salign(overhang=30, gap_penalties_1d=(-450, -50),
           alignment_type='tree', output='ALIGNMENT')

aln.write(file='malign.ali', alignment_format='PIR')
