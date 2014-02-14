from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Gapped

import sys


h1 = open(sys.argv[1], "rU")
h2 = open(sys.argv[2], "rU")

s1 = SeqIO.parse(h1, "pir")
s2 = SeqIO.parse(h2, "pir")

r1 = [r for r in s1]
r2 = [r for r in s2]

align = MultipleSeqAlignment([])
ss = []
#N-Term + C-Term only
for i in r1:
    for j in r2:
        if j.id == i.id:
            d = i.description.split(":")
            d[2] = "FIRST"
            d[4] = "LAST"
            dsc = ":".join(d)
            ss.append(i.id)
            s = i.seq + j.seq
            align.append(SeqRecord(s, id=i.id, name=i.name, description=dsc))
            l = len(s)
            break

#N-Term only
for i in r1:
    if i.id not in ss:
        d = i.description.split(":")
        d[2] = "FIRST "
        d[4] = "LAST "
        dsc = ":".join(d)
        g = l-len(i.seq)
        align.append(SeqRecord(i.seq + "-"*g, id=i.id, name=i.name, description=dsc))

#C-Term only
for j in r2:
    if j.id not in ss:
        d = j.description.split(":")
        d[2] = "FIRST "
        d[4] = "LAST "
        dsc = ":".join(d)
        g = l-len(j.seq)
        align.append(SeqRecord("-"*g + j.seq, id=j.id, name=j.name, description=dsc))


#h = open(sys.argv[3], "w")
#AlignIO.write(align, h, "fasta")
#h.close()
h1.close()
h2.close()


for a in align:
    print ">P1;"+a.id
    print a.description
    l = len(a.seq)
    b = l / 75
    x = l % 75
    
    for i in range(b):
        print a.seq[i*75:i*75+75]

    i += 1
    if x:
        print a.seq[i*75:i*75+75] + "*"

    print ""
