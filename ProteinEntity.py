from Exceptions import *

class ProteinEntity:
    def __init__(self, level, id):
        self.__level = level
        self.__id = id
        self.__full_id = None
        
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
    def add_child(self, child, id):
        if self.has_child_with_id(id):
            raise ProteinManipulationException( \
                "Child %s already exists under %s " \
                % (child, self) )
        elif not child.is_orphaned():
            raise ProteinManipulationException( \
                "Child %s belongs to some other entity %s already " \
                % (child, child.parent()) )
        else:
            self.__children[id] = child
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
        return self.__parent is None

    def set_childrens_parents(self):
        for child in self.__children():
            child.set_parent(self)

    def set_parent(self, parent):
        self.__parent = parent

    def detach_parent(self):
        self.__parent = None


    def full_id(self):
        """Return the full id.

        The full id is a tuple containing all id's starting from
        the top object (Structure) down to the current object. A full id for
        a Residue object e.g. is something like:

        ("1abc", 0, "A", (" ", 10, "A"))

        This corresponds to:

        Structure with id "1abc"
        Model with id 0
        Chain with id "A"
        Residue with id (" ", 10, "A")

        The Residue id indicates that the residue is not a hetero-residue
        (or a water) beacuse it has a blank hetero field, that its sequence
        identifier is 10 and its insertion code "A".
        """
        if self.__full_id==None:
            entity_id=self.id()
            l=[entity_id]
            parent=self.parent()
            while not (parent is None):
                entity_id=parent.id()
                l.append(entity_id)
                parent=parent.parent()
            l.reverse()
            self.__full_id=tuple(l)
        return self.__full_id


class DisorderedProteinEntity(ProteinEntity):
    def __init__(self, level, id):
        ProteinEntity.__init__(self, level, id)
        self.__main_disorder_identifier = None # the atom/residue state with the highest occupancy
        
    def add_child(self, child, id):
        ProteinEntity.add_child(self, child, id)
        child.set_parent(self.parent())
        child.flag_disordered()

    def set_main_disordered_child(self, child):
        self.__main_disorder_identifier = child.disorder_identifier()

    def main_disorder_identifier(self):
        return self.__main_disorder_identifier

    def main_child(self):
        return self.children(self.__main_disorder_identifier)
    # need to check disordered_select on removing child, and to check if disordered_children is empty