"""
Runtime Tests to ensure that everything is working
"""

from Exceptions import *
from PDBInfo import PDBInfo
from ProteinEntity import ProteinEntity
from Structure import Structure
from Model import Model
from Chain import Chain
from Residue import Residue
#from Atom import Atom

info = PDBInfo('1JCL.pdb')
s=Structure('1JCL.pdb')
m=Model(1)
c=Chain(2)
r=Residue(('a','b','c'), 'PHE', 3)
m.add_chain(c)
c.model()
