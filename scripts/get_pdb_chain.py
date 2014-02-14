#This script downloads a PDB file from internet and extract the correspond chain
# Use:
#      python get_pdb_chain.py PDB_ID CHAIN
#
# Note: It will create a direcotry called pdbdb which stores the entire pdb file.
# Carlos Ríos Vera <crosvera@gmail.com> 2014

import sys

from Bio.PDB import *

class ChainSelect(Select):
    def accept_chain(self, chain):
        if chain.id == sys.argv[2]:
            return 1
        else: 
            return 0



if __name__ == "__main__":
    parser = PDBParser()
    
    pdbl = PDBList(pdb="./pdbdb")
    f = pdbl.retrieve_pdb_file(sys.argv[1])

    s = parser.get_structure(sys.argv[1], f)

    io = PDBIO()
    io.set_structure(s)
    io.save(sys.argv[1]+".pdb", ChainSelect())

