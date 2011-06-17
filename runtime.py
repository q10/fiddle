"""
Runtime Tests / API Documentation and Examples
"""

from Structures.Structure import Structure
from Parsers.PerlGenPDBParser import PerlGenPDBParser
from Parsers.RotamerLibParser import RotamerLibParser
from Methods import *
from RotamerLibrary import *
import profile
import datetime

#s=Structure('pdbs/1JCL.pdb')
profile.run("c=PerlGenPDBParser('pdbs/sample.pdb').get_chain()")
r=RotamerLibParser('rotamers/arg.lib')

for i in range(0,232):
    c.remove_residue_with_id(c.residues(i).id())

start_time = datetime.datetime.now()
get_filtered_rotamer_chain(c)
end_time = datetime.datetime.now()
print(end_time - start_time)


"""

# creating new Atom and Residue
new_coords = numpy.array((1.23, 2.56, 3.98), 'f')
a=Atom(12345, 'CG', new_coords)
r=Residue(123, 'PHE')


# print all PDB metadata into a string
s.info()

# hash keys for the PDB metadata hash
s.info_keys()

# access the specific info from hash key
s.info(<key>)

# print info (formatted) to python console
s.show_info()

# get list of all models belonging to structure
s.models()
s.models(<id>)
# <id> format example can bee seen in a model.id()

# get list of all chains belonging to structure
s.chains()

# get list of all residues/disordered residues belonging to structure (no unpacking) (in order of sequence on chain)
s.residues()

# get list of all atoms/disordered atoms belonging to structure (no unpacking) (in order of sequence on residue)
s.atoms()

# add to and remove model from structure
m=Model(1)
s.add_model(m)
s.has_model_with_id(m.id())
s.remove_model_with_id(m.id())

# similar methods, for Models
m = s.models()[0]
m.chains()
m.residues()
m.atoms()

# get chain A from model
m.chains('A')

# methods are 'safe' - if null, remove_chain_with_id will do nothing
c=Chain(2)
m.add_chains(c)
m.has_chain_with_id(c.id())
m.remove_chain_with_id(c.id())

# get the structure this model belongs to
m.structure()

# similar methods, for Chains - residues are in order in relation to chain, and atoms are in order in relation to residue
c=m.chains()[0]
c.residues()
c.atoms()
c.model()
c.structure()

r=Residue(3, 'PHE')
dr=DisorderedResidue(100)
c.add_residue(r)
c.add_residue(dr)
c.has_residue_with_id(r.id())
c.remove_residue_with_id(r.id())


# similar methods, for Residues - atoms are in order in relation to residue
# all these methods work for DisorderedResidues, but the method calls are redirected to the sub-residue with highest occupancy
r = c.residues()[0]
r.atoms()
r.chain()
r.model()
r.structure()

# get the beta carbon of this residue
r.atoms('CB')

r.add_atom(a)
r.has_atom_with_id(a.id())
r.remove_atom_with_id(a.id())

# returns list of atoms, including DisorderedAtoms
r.atoms()

# returns list of atoms, unpacking all atoms from DisorderedAtoms
r.all_atoms()

# returns the residue's sequence starting from 0
# for example, if the residue is the 10th on the chain, sequence_number() will return 9
r.sequence_number()

# residue name - the three-letter name of the residue - ex. phenylalanine => PHE
r.name()

# returns true if it's a disordered residue (occupancy less than 1)
r.is_disordered()


# Methods specific for DisorderedResidues

# returns the residue with highest occupancy
dr.main_residue()

# lists out the residues that occupy this spot
dr.disordered_residues()

# add a new residue in
rr=Residue(4, 'PHE')
dr.add_residue(rr)
dr.has_residue_with_name(rr.id())

# add atom to one of the residues occupying this spot
dr.add_atom_to_residue(a, rr.id())

# remove residue by name() (ex 'PHE'), NOT id()
dr.remove_residue_with_name(rr.name())


# Methods for Atoms, works for DisorderedAtoms as well, but applies to the sub-atom with the highest occupancy

a = r.atoms()[0]
a.chain()
a.model()
a.structure()

# returns only the residue it is under, NOT the DisorderedResidue
a.residue()


# get full numpy array or x y z coordinates
a.coords()
a.x(), a.y(), a.z()

# sets new coords for atom - must be numpy array of length 3
new_coords = numpy.array((1.23, 2.56, 3.98), 'f')
a.set_coords(new_coords)

# get absolute distance from other atom
a.distance_from(other_atom)

# returns element ex. 'C' for carbon
a.element()

# returns true if it's a DisorderedAtom
a.is_disordered()

# the number absolute to its position in a protein
a.serial_number()

# name of atom; ex. 'CG' for the gamma carbon
a.name()

# occupancy of the atom in range (0, 1] -
# ex. a 0.6 means 60% of the time the disordered atom is in the position/state described by this atom.coords()
a.occupancy()


# Methods specific for DisorderedAtoms

da.add_atom(aa)
da.has_atom_with_altloc(aa.altloc())
da.remove_atom_with_altloc(aa.altloc())

# list all the atoms under this DisorderedAtom
da.disordered_atoms()

# returns the Atom under this DisorderedAtom with the highest occupancy
da.main_atom


"""
