"""
Runtime Tests to ensure that everything is working
"""

from Exceptions import *
from StructureBuilder import StructureBuilder
from PDBInfo import PDBInfo
from ProteinEntity import ProteinEntity
from Structure import Structure
from Model import Model
from Chain import Chain
from Residue import *
from Atom import *


s=Structure('1JCL.pdb')
'''
info = PDBInfo('1JCL.pdb')
m=Model(1)
c=Chain(2)
r=Residue(('a','b','c'), 'PHE', 3)
c.add_residue(r)
m.add_chain(c)
s.add_model(m)

rr=Residue(('a','b','c'), 'PHE', 4)
dr=DisorderedResidue(100)
da=DisorderedAtom(101)
c.add_residue(dr)
dr.add_residue(rr)

a=Atom('CA',(1,2,3),'BFACTOR',0.5,'ALTLOC','SERIAL NUM','C')

from Bio.PDB.PDBParser import PDBParser
p=PDBParser(PERMISSIVE=1)
filename="pdb1fat.ent"
s=p.get_structure('1JCL', '1JCL.pdb')
m=s[0]
c=m[0]
r=c[0]
'''
