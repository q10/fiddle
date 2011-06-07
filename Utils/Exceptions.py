"""Enumeration of PDB-based exceptions."""

# General error
class PDBException(Exception):
    pass

# The PDB file cannot be unambiguously represented in the SMCRA
# data structure
class PDBConstructionException(Exception):
    pass

class PDBConstructionWarning(Warning):
    pass

class ProteinManipulationException(Exception):
    pass

class PDBATM2Warning(Warning):
    pass

class GROWSCESCException(Exception):
    pass