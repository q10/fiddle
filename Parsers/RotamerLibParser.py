"""Parser for rotamer library files."""

import os
from itertools import product

class RotamerLibParser:
    """
    Parse a rotamer library file and returns a hash, with residue name capitalized as key and list of tuples as value.

    Each line in a rotamer library is of the following format:
    a da    b db    c dc    ...

    This parser will create a list of lists
    [[a, a+da, a-da], [b, b+db, b-db], [c, c+dc, d-dc], ...]

    The product method will generate the cartesian product of all the sublists as tuples,
    which will be inserted into the rotamer library list, self.__all_torsions
    """

    def __init__(self, filename):
        assert(isinstance(filename, str))
        self.__residue_name = os.path.basename(filename).split('.')[0].upper()
        self.__num_columns = 0
        self.__all_torsions = []
        self._parse_rotamers(open(filename).readlines())

    def residue_name(self):
        return str(self.__residue_name)

    def torsions(self):
        return self.__all_torsions

    def _parse_rotamers(self, lines):
        self.__num_columns = len(lines[0].split())
        assert(self.__num_columns % 2 == 0)
        self.__num_columns = int(self.__num_columns / 2)

        for i in range(len(lines)):
            torsions = []
            items = lines[i].split()
            for k in range(0, len(items), 2):
                a = float(items[k])
                b = float(items[k+1])
                torsions.append([a, a+b, a-b])
            assert(len(torsions) == self.__num_columns)

            for cartesian_product in product(*torsions):
                assert(len(cartesian_product) == self.__num_columns)
                self.__all_torsions.append(cartesian_product)