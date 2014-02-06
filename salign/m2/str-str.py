
from modeller import *
import sys

log.verbose
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

structures = (('2vqp', 'A'),                               #RSV
              ('4lp7', 'A'),                               #Metapneumovirus
              ('4ldd', 'B'), ('1es6', 'A'), ('4ldi', 'A'), #Ebola
              ('4ld8', 'A'), ('3tcq', 'A'), ('4ldb', 'C'), #Ebola
              ('4g1l', 'B'), ('4g1g', 'B'), ('4g1o', 'B')) #Newcastle

#for (_code, _start, _end, _code_ap) in (('1BL0', 'FIRST:A', 'LAST:A', 'A'), ('1XS9', 'FIRST:A', 'LAST:A', 'A'), ('1D5Y', 'FIRST:C', 'LAST:C', 'C'), ('3OOU', 'FIRST:A', 'LAST:A', 'A'), ('3LSG', 'FIRST:A', 'LAST:A', 'A'), ('3OIO', 'FIRST:A', 'LAST:A', 'A'), ('3MKL', 'FIRST:B', 'LAST:B', 'B'), ('3MN2', 'FIRST:B', 'LAST:B', 'B'), ('2K9S', 'FIRST:A', 'LAST:A', 'A'), ('3GBG', 'FIRST:A', 'LAST:A', 'A')):
#   mdl = model(env, file=_code, model_segment=(_start, _end))
#   aln.append_model(mdl, atom_files=_code, align_codes=_code+_code_ap)
for (code, chain) in structures:
    mdl = model(env, file=code, model_segment=('FIRST:'+chain, 'LAST:'+chain))
    aln.append_model(mdl, atom_files=code, align_codes=code+chain)

aln.write(file='str-str_out.aliIni', alignment_format='pir')
fil = "str-str_out.aliIni"

opfile = "str-str_out.aliMid"
opfile2 = "str-str_out.ali"

#opfile = "salign_local_mid.ali"
#opfile1 = "salign_local.pap"
#opfile2 = "salign_local.ali"

nejon = True
poi = False

#log.verbose
#env = environ()
#env.io.atom_files_directory = 'structures/'


def frange(start,end=None,inc=None):
#  "A range function that accepts floating point increments"

  if end == None:
    end = start + 0.0
    start = 0.0
  else:
    start +=0.0

  if inc == None:
    inc = 1.0

  count = int((end - start)/inc)
  if start + (count*inc) != end:
    count += 1


  L = [None,]*count
  for i in xrange(count):
    L[i] = start + i*inc

  return L


# -- Script that takes in user specified feature weights and gap_penalties_1d, 
# -- given an input alignment
def salign_fw_gaps1(aln,fil,fw,ogp,egp):

 log.verbose
 env = environ()
 env.io.atom_files_directory = 'structures/'
# aln = alignment(env)
# aln.append(file=fil, align_codes='all')
 nseg = 2


 L =  aln.salign(rms_cutoff=3.5,
      normalize_pp_scores=False,
      rr_file='$(LIB)/as1.sim.mat', overhang=0,
      auto_overhang=True, overhang_auto_limit=5, overhang_factor=1,
      gap_penalties_1d=(ogp, egp),
#     local_alignment=True, matrix_offset = -0.2,
      local_alignment=False, matrix_offset = -0.2,
      gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
#     write_weights=False, output_weights_file ='salign.wgt',
      dendrogram_file='str-str.tree',
      alignment_type='tree',
      nsegm=nseg,
      feature_weights=fw,
      improve_alignment=True, fit=True, write_fit=False , write_whole_pdb=False,
      output='ALIGNMENT QUALITY' )

 return L

# -- Script that takes in user specified feature weights and gap_penalties_3d, 
# -- given an input alignment
def salign_fw_gaps3(aln,fil,fw,ogp3d,egp3d,wf):

 log.verbose
 env = environ()
 env.io.atom_files_directory = 'structures/'
# aln = alignment(env)
# aln.append(file=fil, align_codes='all')
 nseg = 2
 ogp = ogp3d
 egp = egp3d


 L = aln.salign(rms_cutoff=3.5,
      normalize_pp_scores=False,
      rr_file='$(LIB)/as1.sim.mat', overhang=0,
      auto_overhang=True, overhang_auto_limit=5, overhang_factor=1,
      gap_penalties_1d=(ogp, egp),
#     local_alignment=True, matrix_offset = -0.2,
      local_alignment=False, matrix_offset = -0.2,
      gap_penalties_3d=(ogp3d, egp3d), gap_gap_score=0, gap_residue_score=0,
#     write_weights=False, output_weights_file ='salign.wgt',
      dendrogram_file='str-str.tree',
      alignment_type='tree',
      nsegm=nseg,
      feature_weights=fw,
#      improve_alignment=True, fit=True, write_fit=wf,  write_whole_pdb=False,
      improve_alignment=True, fit=True, write_fit=True,  write_whole_pdb=False,
      output='ALIGNMENT QUALITY' )

 return L

# -- Script that takes in user specified feature weights and gap_penalties_1d, 
# -- given an input alignment
def salign_fw_local_gaps1(aln,fil,fw,ogp,egp,mat_off):

 log.verbose
 env = environ()
 env.io.atom_files_directory = 'structures/'
# aln = alignment(env)
# aln.append(file=fil, align_codes='all')
 nseg = 2


 L =  aln.salign(rms_cutoff=3.5,
      normalize_pp_scores=False,
      rr_file='$(LIB)/as1.sim.mat', overhang=0,
#     auto_overhang=True, overhang_auto_limit=5, overhang_factor=1,
      gap_penalties_1d=(ogp, egp),
      local_alignment=True, matrix_offset = mat_off, matrix_offset_3d = -0.5,
#     local_alignment=False, matrix_offset = -0.2,
      gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
#     write_weights=False, output_weights_file ='salign.wgt',
      dendrogram_file='str-str.tree',
      alignment_type='tree',
      nsegm=nseg,
      feature_weights=fw,
      improve_alignment=True, fit=True, write_fit=False ,
      output='ALIGNMENT QUALITY' )

 return L

# -- Script that takes in user specified feature weights and gap_penalties_3d, 
# -- given an input alignment
def salign_fw_local_gaps3(aln,fil,fw,ogp3d,egp3d,mat_off,mat_off_3d,wf):

 log.verbose
 env = environ()
 env.io.atom_files_directory = 'structures/'
# aln = alignment(env)
# aln.append(file=fil, align_codes='all')
 nseg = 2
 ogp = ogp3d
 egp = egp3d


 L = aln.salign(rms_cutoff=3.5,
      normalize_pp_scores=False,
      rr_file='$(LIB)/as1.sim.mat', overhang=0,
#     auto_overhang=True, overhang_auto_limit=5, overhang_factor=1,
      gap_penalties_1d=(ogp, egp),
      local_alignment=True, matrix_offset = mat_off, matrix_offset_3d = mat_off_3d,
#     local_alignment=False, matrix_offset = -0.2,
      gap_penalties_3d=(ogp3d, egp3d), gap_gap_score=0, gap_residue_score=0,
#     write_weights=False, output_weights_file ='salign.wgt',
      dendrogram_file='str-str.tree',
      alignment_type='tree',
      nsegm=nseg,
      feature_weights=fw,
#      improve_alignment=True, fit=True, write_fit=wf,
      improve_alignment=True, fit=True, write_fit=True,
      output='ALIGNMENT QUALITY' )

 return L


# -- Iterating over values of gap penalties and nsegm
qmax = 0.0
nsegm = 2
fw1=(1., 0., 0., 0., 1., 0.)
fw2=(0., 1., 0., 0., 0., 0.)
fw3=(0., 0., 0., 0., 1., 0.)

# -- Iterating over gap penalties 1D to get initial alignments
for ogp in frange(-150,1,30):
   for egp in frange(-50,1,10):
      for mo in frange(-3.0, -0.05, 0.3) :
          aln = alignment(env)
          aln.append(file=fil, align_codes='all')
          try:
            qwlty1 = salign_fw_local_gaps1(aln,fil,fw1,ogp,egp,mo)
            if qwlty1.qscorepct >= qmax:
   	 	qmax = qwlty1.qscorepct
                aln.write(file=opfile, alignment_format='PIR')
		win_ogp = ogp
		win_egp = egp
		win_mo = mo
            print "Qlty scrs", ogp,"\t",egp,"\t",qwlty1.qscorepct
          except ModellerError, detail:
            print "Set of parameters",fw1,ogp,egp,"resulted in the following error\t"+str(detail)
          del(aln)


# -- Iterating over gap panelties 3D to get final alignments
for ogp3d in frange(0,3,1) :
   for egp3d in range (2,5,1) :
            aln = alignment(env)
            aln.append(file=opfile, align_codes='all')
            try:
               qwlty2 = salign_fw_gaps3(aln,opfile,fw2,ogp3d,egp3d,poi)
	       if qwlty2.qscorepct >= qmax:
                  qmax = qwlty2.qscorepct
#                  aln.write(file=opfile1, alignment_format='PAP')
                  aln.write(file=opfile2, alignment_format='PIR')
		  win_ogp3d = ogp3d
		  win_egp3d = egp3d
               print "Qlty scrs", ogp3d,"\t",egp3d,"\t",qwlty2.qscorepct
            except ModellerError,detail:
               print "Set of parameters",fw2,ogp3d,egp3d,"resulted in the following error\t"+str(detail)
            del (aln)

#print "final max quality = ",qmax

# -- try alternate initial alignments only if the qmax score is less than 70%

# - ******** qmax threshold for additional iterations to be determined *******
qmax_old = qmax
if (qmax_old <= 70) :
#  qmax = 0.0
   for ogp in frange(0.0,2.2,0.3):
      for egp in frange(0.1,2.3,0.3):
         for mo in frange (-3.0, -0.05, 0.3) :
             aln = alignment(env)
             aln.append(file=fil, align_codes='all')
             try:
               qwlty1 = salign_fw_local_gaps1(aln,fil,fw3,ogp,egp,mo)
               if qwlty1.qscorepct >= qmax:
                    qmax = qwlty1.qscorepct
                    aln.write(file=opfile, alignment_format='PIR')
		    win_ogp = ogp
		    win_egp = egp
		    win_mo = mo
               print "Qlty scrs", ogp,"\t",egp,"\t",qwlty1.qscorepct
             except ModellerError, detail:
               print "Set of parameters",fw3,ogp,egp,"resulted in the following error\t"+str(detail)
             del(aln)

# -- Iterating over gap panelties 3D to get final alignments
   for ogp3d in frange(0,3,1) :
      for egp3d in range (2,5,1) :

               aln = alignment(env)
               aln.append(file=opfile, align_codes='all')
               try:
                  qwlty2 = salign_fw_gaps3(aln,opfile,fw2,ogp3d,egp3d,poi)
                  if qwlty2.qscorepct >= qmax:
                     qmax = qwlty2.qscorepct
#                     aln.write(file=opfile1, alignment_format='PAP')
                     aln.write(file=opfile2, alignment_format='PIR')
		     win_ogp3d = ogp3d
		     win_egp3d = egp3d
                  print "Qlty scrs", ogp3d,"\t",egp3d,"\t",qwlty2.qscorepct
               except ModellerError,detail:
                  print "Set of parameters",fw2,ogp3d,egp3d,"resulted in the following error\t"+str(detail)
               del (aln)

print "final max quality = ",qmax


aln = alignment(env)
aln.append(file=opfile, align_codes = 'all')
salign_fw_gaps3(aln,opfile,fw2,win_ogp3d,win_egp3d,nejon)
print "Completed successfully"


# Make two comparisons: no cutoffs, and 3.5A/60 degree cutoffs for RMS, DRMS,
# and dihedral angle comparisons:
aln.compare_structures(rms_cutoffs=[999]*11)
aln.compare_structures(rms_cutoffs=(3.5, 3.5, 60, 60, 60, 60, 60, 60, 60, 60, 60))
print "Complete compare structures"
