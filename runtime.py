"""
Runtime Tests to ensure that everything is working
"""

from Exceptions import *
from ProteinEntity import ProteinEntity
from Structure import Structure
from Model import Model
from Chain import Chain
#from Residue import Residue
#from Atom import Atom

c=Chain(2)
m=Model(1)
m.add_chain(c)
c.model()
