# Illustrates the SALIGN iteractive multiple structure alignment
from modeller import *
import modeller.salign

log.none()
env = environ()
env.io.atom_files_directory = ['.', '../atom_files']

aln = alignment(env)
structures = (('4lp7', 'A'), ('4hi6', 'C'), ('4ldd', 'B'), 
              ('1h2d', 'A'), ('1h2c', 'A'), ('1es6', 'A'), 
              ('4ldi', 'A'), ('4ldm', 'A'), ('4hiu', 'A'), 
              ('3fij', 'A'), ('4hi5', 'A'), ('4hiw', 'A'), 
              ('4hiy', 'A'), ('4ld8', 'A'), ('4hit', 'C'), 
              ('3tcq', 'A'), ('4ldb', 'C'), ('4g1l', 'B'), 
              ('4g1g', 'B'), ('4g1o', 'B'))

for (code, chain) in structures:
    mdl = model(env, file=code, model_segment=('FIRST:'+chain, 'LAST:'+chain))
    aln.append_model(mdl, atom_files=code, align_codes=code+chain)

modeller.salign.iterative_structural_align(aln)

aln.write(file='2vqp-like.pap', alignment_format='PAP')
aln.write(file='2vqp-like.ali', alignment_format='PIR')
