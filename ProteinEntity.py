class ProteinEntity:
    def __init__(self):
        self.__parent = None

    def level(self):
        return self.level

    def set_parent(self, parent):
        self.__parent = parent
