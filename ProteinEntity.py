from Exceptions import *

class ProteinEntity:
    def __init__(self, level, id):
        self.__level = level
        self.__id = id
        
        self.__parent = None
        self.__children = {}

    def children(self, hash_key=None):
        if hash_key is None:
            return list(self.__children.values())
        else:
            try:
                return self.__children[hash_key]
            except:
                return None

    def id(self):
        return self.__id

    def level(self):
        return self.__level



    # Sub-entity manipulation methods
    def add_child(self, child):
        if self.has_child_with_id(child.id()):
            raise ProteinManipulationException( \
                "Child %s already exists under %s " \
                % (child, self) )
        elif child.is_orphaned():
            raise ProteinManipulationException( \
                "Child %s belongs to some other entity %s already " \
                % (child, child.parent()) )
        else:
            self.__children[child.id()] = child
            child.set_parent(self)

    def has_child_with_id(self, id):
        try:
            child = self.__children[id]
        except:
            return False
        else:
            return True

    def remove_child_with_id(self, id):
        if self.has_child_with_id(id):
            self.__children[id].detach_parent()
            del self.__children[id]
        else:
            raise ProteinManipulationException( \
                "Child with id %s does not exist under %s " \
                % (id, self) )

    def parent(self):
        return self.__parent
    
    def is_orphaned(self):
        return self.__parent is not None

    def set_childrens_parents(self):
        for child in self.__children():
            child.set_parent(self)

    def set_parent(self, parent):
        self.__parent = parent

    def detach_parent(self):
        self.__parent = None

