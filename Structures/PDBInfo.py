from Parsers.PDBHeaderParser import parse_pdb_header

class PDBInfo:
    def __init__(self, filename):
        file = open(filename,'r')
        self.__info_table = parse_pdb_header(file)
        file.close()

    def __str__(self):
        return self.full_description()

    def keys(self):
        return list(self.__info_table.keys())

    def full_description(self):
        full_description = ""
        for k, v in self.__info_table.items():
            full_description += ("-"*40) + '\n'
            full_description += str(k) + '\n'
            full_description += str(v) + '\n'
        return full_description

    def __getitem__(self, hash_key):
        try:
            return self.__info_table[hash_key]
        except:
            return None