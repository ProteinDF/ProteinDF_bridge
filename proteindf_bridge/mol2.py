#!/usr/bin/env python

import logging
logger = logging.getLogger(__name__)

from .position import Position
from .atom import Atom
from .atomgroup import AtomGroup
from .select import Select_Path

class SimpleMol2(object):
    """
    """
    def __init__(self, atomgroup = AtomGroup()):
        self._atomgroup = AtomGroup()
        self._atom_index_table = []
        self.set_by_atomgroup(atomgroup)


    def set_by_atomgroup(self, atomgroup):
        assert(isinstance(atomgroup, AtomGroup))
        self._atomgroup = atomgroup

        # make atom index table
        self._atom_index_table = [None] * self._atomgroup.get_number_of_all_atoms()
        index = 0
        for atom in self._atomgroup.get_atom_list():
            self._atom_index_table[index] = atom.name
            index += 1


    def save(self, file_path):
        with open(file_path, "w") as f:
            contents = str(self)
            f.write(contents)


    def _get_contents_molecule(self):
        """molecule section

        mol_name
        num_atoms [num_bonds [num_subst [num_feat [num_sets]]]]
        mol_type 
        charge_type
        [status_bits 
        [mol_comment]]

        mol_type (string) = the molecule type: SMALL, BIOPOLYMER, PROTEIN, NUCLEIC_ACID, SACCHARIDE
        charge_type (string) = the type of charges associated with the molecule:
            NO_CHARGES, DEL_RE, GASTEIGER, GAST_HUCK, HUCKEL, PULLMAN, GAUSS80_CHARGES, AMPAC_CHARGES, MULLIKEN_CHARGES, DICT_ CHARGES, MMFF94_CHARGES, USER_CHARGES
        """
        molecular_name = self._atomgroup.name
        num_atoms = self._atomgroup.get_number_of_all_atoms()
        mol_type = "SMALL"
        charge_type = "USER_CHARGES"

        answer = ""
        answer += "@<TRIPOS>MOLECULE\n"
        answer += "{}\n".format(molecular_name)
        answer += "{num_atoms}\n".format(num_atoms=num_atoms)
        answer += "{mol_type}\n".format(mol_type=mol_type)
        answer += "{charge_type}\n".format(charge_type=charge_type)
        answer += "\n"

        return answer

    def _get_contents_atom(self):
        """
        format: 
        atom_id atom_name x y z atom_type [subst_id [subst_name [charge [status_bit]]]]
        """

        answer = ""
        answer += "@<TRIPOS>ATOM\n"
        atom_id = 1
        for atom_key, atom in self._atomgroup.atoms():
            #print(atom)
            atom_name = atom.name
            x = atom.position.x
            y = atom.position.y
            z = atom.position.z
            atom_type = atom.symbol
            answer += "{atom_id:<5} {atom_name:<3} {x: 8.3f} {y: 8.3f} {z: 8.3f} {atom_type}\n".format(
                atom_id=atom_id, atom_name=atom_name,
                x=x, y=y, z=z, atom_type=atom_type
                )
            atom_id += 1
        answer += "\n"

        return answer

    def _get_contents_bond(self):
        """
        format: 
        bond_id origin_atom_id target_atom_id bond_type [status_bits]
        """
        answer = ""
        answer += "@<TRIPOS>BOND\n"

        bond_id = 1
        for bond_info in self._atomgroup.get_bond_list():
            (atom_path1, atom_path2, bond_order) = bond_info
            
            selector_atom1 = Select_Path(atom_path1)
            selector_atom2 = Select_Path(atom_path2)
            ag1 = self._atomgroup.select(selector_atom1)
            ag2 = self._atomgroup.select(selector_atom2)
            atom1 = ag1.get_atom_list()[0]
            atom2 = ag2.get_atom_list()[0]
            #print(atom1.name, atom2.name)

            atom_id1 = self._atom_index_table.index(atom1.name)
            atom_id2 = self._atom_index_table.index(atom2.name)
            bond_type = bond_order
            answer += "{bond_id:<5} {atom_id1:<5} {atom_id2:<5} {bond_type}\n".format(
                bond_id=bond_id, atom_id1=atom_id1, atom_id2=atom_id2,
                bond_type=bond_type
                )
            bond_id += 1
        answer += "\n"

        return answer
    
    def __str__(self):
        answer = ""
        answer += self._get_contents_molecule()
        answer += self._get_contents_atom()
        answer += self._get_contents_bond()
        return answer
