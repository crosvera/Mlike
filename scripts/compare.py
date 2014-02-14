
from modeller import *
import sys

log.verbose()
env = environ()
env.io.atom_files_directory = '../../pdbs/'
aln = alignment(env)
#structures = (('4lp7', 'A'), ('4hi6', 'C'), ('4ldd', 'B'), 
#              ('1h2d', 'A'), ('1h2c', 'A'), ('1es6', 'A'), 
#              ('4ldi', 'A'), ('4ldm', 'A'), ('4hiu', 'A'), 
#              ('3f1j', 'A'), ('4hi5', 'A'), ('4hiw', 'A'), 
#              ('4hiy', 'A'), ('4ld8', 'A'), ('4hit', 'C'), 
#              ('3tcq', 'A'), ('4ldb', 'C'), ('4g1l', 'B'), 
#              ('4g1g', 'B'), ('4g1o', 'B'), ('2ykd', 'A'))

#structures = (('2vqp', '129:A', '259:A', 'A'),      #RSV
#              ('4lp7', '138:A', '254:A', 'A'),                              #Metapneumovirus
#              ('4ldd', '191:B', '321:B', 'B'), ('1es6', '191:A', '321:A', 'A'), ('4ldi', '191:A', '321:A', 'A'), #Ebola
#              ('4ld8', '191:A', '321:A', 'A'), ('3tcq', '191:A', '321:A', 'A'), ('4ldb', '191:C', '321:C', 'C'), #Ebola
#              ('4g1l', '182:B', '252:B', 'B'), ('4g1g', '182:B', '252:B', 'B'), ('4g1o', '182:B', '252:B', 'B')) #Newcastle

#structures = (('2vqp', '127:A', 'LAST:A', 'A'), ('2ykd', '131:A', 'LAST:A', 'A'),     #RSV
#              ('4lp7', '137:D', 'LAST:D', 'D'),                              #Metapneumovirus
#              ('4ldd', '191:B', 'LAST:B', 'B'), ('1es6', '191:A', 'LAST:A', 'A'), ('4ldi', '191:A', 'LAST:A', 'A'), #Ebola
#              ('4ld8', '191:A', 'LAST:A', 'A'), ('3tcq', '191:A', 'LAST:A', 'A'), ('4ldb', '191:C', 'LAST:C', 'C'), #Ebola
#              ('4g1l', '189:B', 'LAST:B', 'B'), ('4g1g', '189:B', 'LAST:B', 'B'), ('4g1o', '189:B', 'LAST:B', 'B'), #Newcastle
#              ('4hi6', 'FIRST:C', 'LAST:C', 'C'), ('4hiu', 'FIRST:A', 'LAST:A', 'A'), ('3f1j', 'FIRST:A', 'LAST:A', 'A'), #BDV
#              ('4hi5', 'FIRST:A', 'LAST:A', 'A'), ('4hiw', 'FIRST:A', 'LAST:A', 'A'), ('4hiy', 'FIRST:A', 'LAST:A', 'A'), #BDV
#              ('4hit', 'FIRST:C', 'LAST:C', 'C')) #BDV


#for (_code, _start, _end, _code_ap) in (('1BL0', 'FIRST:A', 'LAST:A', 'A'), ('1XS9', 'FIRST:A', 'LAST:A', 'A'), ('1D5Y', 'FIRST:C', 'LAST:C', 'C'), ('3OOU', 'FIRST:A', 'LAST:A', 'A'), ('3LSG', 'FIRST:A', 'LAST:A', 'A'), ('3OIO', 'FIRST:A', 'LAST:A', 'A'), ('3MKL', 'FIRST:B', 'LAST:B', 'B'), ('3MN2', 'FIRST:B', 'LAST:B', 'B'), ('2K9S', 'FIRST:A', 'LAST:A', 'A'), ('3GBG', 'FIRST:A', 'LAST:A', 'A')):
#for (_code, _start, _end, _code_ap) in structures:
#   mdl = model(env, file=_code, model_segment=(_start, _end))
#   aln.append_model(mdl, atom_files=_code, align_codes=_code+_code_ap)

aln.append(file=sys.argv[1], align_codes='all')
aln.compare_structures(fit=True, rms_cutoffs=[999]*11)
aln.compare_structures(fit=True)
print "Complete compare structures"
aln.id_table('id.mat')
