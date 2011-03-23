class ProteinEntity:
    def __init__(self):
        self.__parent = None
        self.__children = {}

    def level(self):
        return self.level

    def __children(self, hash_key=None):
        if index is None:
            return list(self.__children.values())
        else:
            return self.__children[hash_key]

    def set_childrens_parents(self):
        for child in self.__children():
            child.set_parent(self)
            
    def set_parent(self, parent):
        self.__parent = parent
