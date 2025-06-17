from typing import Any, Callable, Dict, List, Optional, Tuple, Iterable

from .store import Store
from .vector3d import Vector3d, pos_distance

radius_by_elem = {
    "H": 1.20,
    "N": 1.55,
    "NA": 2.27,
    "CU": 1.40,
    "CL": 1.75,
    "C": 1.70,
    "O": 1.52,
    "I": 1.98,
    "P": 1.80,
    "B": 1.85,
    "BR": 1.85,
    "S": 1.80,
    "SE": 1.90,
    "F": 1.47,
    "FE": 1.80,
    "K": 2.75,
    "MN": 1.73,
    "MG": 1.73,
    "ZN": 1.39,
    "HG": 1.8,
    "XE": 1.8,
    "AU": 1.8,
    "LI": 1.8,
    ".": 1.8,
}


# Utility functions


def push_to_list_in_dict(dict_obj: Dict, key: Any, value: Any):
    """Add value to list in dictionary, creating list if key doesn't exist"""
    if key not in dict_obj:
        dict_obj[key] = []
    dict_obj[key].append(value)


def get_index_from_values(
    values: List, value: Any, is_equal_fn: Optional[Callable] = None
) -> int:
    """Get index of value in table, adding if not present"""
    if is_equal_fn is None:
        if value not in values:
            values.append(value)
        return values.index(value)
    else:
        for i, this_val in enumerate(values):
            if is_equal_fn(value, this_val):
                return i
        values.append(value)
        return len(values) - 1


def int_to_bool(i: int) -> bool:
    """Convert integer to boolean"""
    return i == 1


def bool_to_int(b: bool) -> int:
    """Convert boolean to integer"""
    return 1 if b else 0


def int_to_char(i: int) -> str:
    """Convert integer to character"""
    return chr(i) if i else ""


def char_to_int(c: str) -> int:
    """Convert character to integer"""
    return ord(c) if c else 0


def in_array(item: Any, array: List) -> bool:
    """Check if item is in array"""
    return item in array


# Store field definitions
ATOM_STORE_FIELDS = [
    ("x", "float32"),
    ("y", "float32"),
    ("z", "float32"),
    ("bfactor", "float32"),
    ("occupancy", "float32"),
    ("alt", "uint8"),
    ("i_atom_type", "uint16"),
    ("i_elem", "uint16"),
    ("i_res", "uint32"),
    ("i_chain", "int32"),
    ("radius", "float32"),
    ("model", "int8"),
]

RESIDUE_STORE_FIELDS = [
    ("atom_offset", "uint32"),
    ("atom_count", "uint16"),
    ("i_central_atom", "uint32"),
    ("i_res_type", "uint16"),
    ("i_chain", "uint8"),
    ("i_structure", "int8"),
    ("res_num", "int32"),
    ("ins_code", "uint8"),
    ("ss", "uint8"),
    ("i_color", "uint16"),
    ("is_polymer", "uint8"),
]


class AtomProxy:
    """Proxy for accessing atom data from stores"""

    def __init__(self, soup: "Soup", i_atom: Optional[int] = None):
        self.soup = soup
        self.i_atom = None
        self._pos = Vector3d()
        if isinstance(i_atom, int):
            self.load(i_atom)

    def load(self, i_atom: int):
        """Load atom at given index"""
        self.i_atom = i_atom
        self._pos.x = self.soup.atom_store.x[self.i_atom]
        self._pos.y = self.soup.atom_store.y[self.i_atom]
        self._pos.z = self.soup.atom_store.z[self.i_atom]
        return self

    @property
    def pos(self) -> Vector3d:
        """Get atom position"""
        return self._pos

    @property
    def bfactor(self) -> float:
        """Get B-factor"""
        return self.soup.atom_store.bfactor[self.i_atom]

    @bfactor.setter
    def bfactor(self, bfactor: float):
        """Set residue index"""
        self.soup.atom_store.bfactor[self.i_atom] = bfactor

    @property
    def occupancy(self) -> float:
        """Get B-factor"""
        return self.soup.atom_store.occupancy[self.i_atom]

    @occupancy.setter
    def occupancy(self, occupancy: float):
        """Set residue index"""
        self.soup.atom_store.occupancy[self.i_atom] = occupancy

    @property
    def elem(self) -> str:
        """Get element"""
        i_elem = self.soup.atom_store.i_elem[self.i_atom]
        return self.soup.elems[i_elem]

    @property
    def atom_type(self) -> str:
        """Get atom type"""
        i_atom_type = self.soup.atom_store.i_atom_type[self.i_atom]
        return self.soup.atom_types[i_atom_type]

    @property
    def alt(self) -> str:
        """Get alternate location"""
        return int_to_char(self.soup.atom_store.alt[self.i_atom])

    @alt.setter
    def alt(self, c: str):
        """Set alternate location"""
        self.soup.atom_store.alt[self.i_atom] = char_to_int(c)

    @property
    def i_res(self) -> int:
        """Get residue index"""
        return self.soup.atom_store.i_res[self.i_atom]

    @i_res.setter
    def i_res(self, i_res: int):
        """Set residue index"""
        self.soup.atom_store.i_res[self.i_atom] = i_res

    @property
    def radius(self) -> float:
        """Get atom radius"""
        return self.soup.atom_store.radius[self.i_atom]

    @property
    def model(self) -> int:
        """Get model number"""
        return self.soup.atom_store.model[self.i_atom]

    @model.setter
    def model(self, model: int):
        """Set model number"""
        self.soup.atom_store.model[self.i_atom] = model

    def __repr__(self):
        """String representation of atom"""
        if self.i_atom is None:
            return "AtomProxy(unloaded)"

        res_proxy = self.soup.get_residue_proxy(self.i_res)
        return (
            f"{res_proxy.chain}{res_proxy.res_num}{res_proxy.ins_code.strip()}-{res_proxy.res_type}:{self.atom_type}"
        )



class ResidueProxy:
    """Proxy for accessing residue data from stores"""

    def __init__(self, soup: "Soup", i_res: Optional[int] = None):
        self.soup = soup
        self.i_res = None
        if isinstance(i_res, int):
            self.load(i_res)

    def load(self, i_res: int):
        """Load residue at given index"""
        self.i_res = i_res
        return self

    @property
    def i_atom(self) -> int:
        """Get central atom index"""
        return self.soup.residue_store.i_central_atom[self.i_res]

    @i_atom.setter
    def i_atom(self, i_atom: int):
        """Set central atom index"""
        self.soup.residue_store.i_central_atom[self.i_res] = i_atom

    @property
    def i_chain(self) -> int:
        """Get chain index"""
        return self.soup.residue_store.i_chain[self.i_res]

    @property
    def chain(self) -> str:
        """Get chain identifier"""
        return self.soup.chains[self.i_chain]

    @property
    def i_structure(self) -> int:
        """Get structure index"""
        return self.soup.residue_store.i_structure[self.i_res]

    @i_structure.setter
    def i_structure(self, i_structure: int):
        """Set structure index"""
        self.soup.residue_store.i_structure[self.i_res] = i_structure

    @property
    def res_num(self) -> int:
        """Get residue number"""
        return self.soup.residue_store.res_num[self.i_res]

    @res_num.setter
    def res_num(self, res_num: int):
        """Set insertion code"""
        self.soup.residue_store.res_num[self.i_res] = res_num

    @property
    def ins_code(self) -> str:
        """Get insertion code"""
        return int_to_char(self.soup.residue_store.ins_code[self.i_res])

    @ins_code.setter
    def ins_code(self, c: str):
        """Set insertion code"""
        self.soup.residue_store.ins_code[self.i_res] = char_to_int(c)

    @property
    def is_polymer(self) -> bool:
        """Check if residue is polymer"""
        return int_to_bool(self.soup.residue_store.is_polymer[self.i_res])

    @is_polymer.setter
    def is_polymer(self, v: bool):
        """Set polymer flag"""
        self.soup.residue_store.is_polymer[self.i_res] = bool_to_int(v)

    @property
    def res_type(self) -> str:
        """Get residue type"""
        i_res_type = self.soup.residue_store.i_res_type[self.i_res]
        return self.soup.res_types[i_res_type]

    @property
    def ss(self) -> str:
        """Get secondary structure"""
        return int_to_char(self.soup.residue_store.ss[self.i_res])

    @ss.setter
    def ss(self, c: str):
        """Set secondary structure"""
        self.soup.residue_store.ss[self.i_res] = char_to_int(c)

    def get_atom_indices(self) -> List[int]:
        """Get indices of atoms in this residue"""
        i_start = self.soup.residue_store.atom_offset[self.i_res]
        n = self.soup.residue_store.atom_count[self.i_res]
        i_end = i_start + n
        return list(range(i_start, i_end))

    def get_i_atom(self, atom_type: str) -> Optional[int]:
        """Get atom index for given atom type"""
        for i_atom in self.get_atom_indices():
            i_atom_type = self.soup.atom_store.i_atom_type[i_atom]
            test_atom_type = self.soup.atom_types[i_atom_type]
            if test_atom_type == atom_type:
                return i_atom
        return None

    def check_atom_types(self, atom_types: List[str]) -> bool:
        """Check if residue has all specified atom types"""
        for atom_type in atom_types:
            if self.get_i_atom(atom_type) is None:
                return False
        return True


class Soup:
    """
    Soup: main data object that holds information
    about protein structure. The soup will be
    embedded in a SoupView that will handle
    all the different viewing options.
    Allowable mutations on the Soup
    will be made via the Controller.
    """

    def __init__(self):
        self.parsing_error = ""
        self.title = ""
        self.structure_ids = []
        self.structure_id = None
        self.i_structure = -1
        self.chains = []

        # Data stores
        self.atom_store = Store(ATOM_STORE_FIELDS)
        self.residue_store = Store(RESIDUE_STORE_FIELDS)

        # Lookup tables
        self.res_ids = []
        self.residue_normal = {}
        self.count_res_added = 0

        # Proxies for data access
        self.residue_proxy = ResidueProxy(self)
        self.atom_proxy = AtomProxy(self)

        # Value tables
        self.elems = []
        self.atom_types = []
        self.res_types = []

        self.max_length = None

    def is_empty(self) -> bool:
        """Check if soup has no atoms"""
        return self.get_atom_count() == 0

    def push_structure_id(self, structure_id: str, title: str = None):
        """Add structure ID"""
        self.structure_id = structure_id
        self.structure_ids.append(structure_id)
        self.i_structure = self.structure_ids.index(structure_id)

        if title:
            id_upper = self.structure_id.upper()
            self.title = f"""[<a href="http://www.rcsb.org/structure/{id_upper}">{id_upper}</a>]{title}"""

    def calc_atom_configuration(self):
        self.max_length = self.calc_max_length()

    def add_atom(
        self,
        x: float,
        y: float,
        z: float,
        bfactor: float,
        alt: str,
        atom_type: str,
        elem: str,
        res_type: str,
        res_num: int,
        ins_code: str,
        chain: str,
        occupancy: float = 1.0,
        model: int = 1,
    ):
        i_atom = self.atom_store.count
        self.atom_store.increment()

        self.atom_store.x[i_atom] = x
        self.atom_store.y[i_atom] = y
        self.atom_store.z[i_atom] = z

        self.atom_store.bfactor[i_atom] = bfactor
        self.atom_store.occupancy[i_atom] = occupancy
        self.atom_store.alt[i_atom] = char_to_int(alt)

        self.atom_store.i_atom_type[i_atom] = get_index_from_values(
            self.atom_types, atom_type
        )

        self.atom_store.i_elem[i_atom] = get_index_from_values(self.elems, elem)
        self.atom_store.radius[i_atom] = radius_by_elem.get(elem, radius_by_elem["."])

        n_res = self.get_residue_count()

        is_new_res = False
        if n_res == 0:
            is_new_res = True
        else:
            # This would use actual residue proxy in real implementation
            last_res_idx = n_res - 1
            if (
                self.residue_store.res_num[last_res_idx] != res_num
                or self.residue_store.ins_code[last_res_idx] != char_to_int(ins_code)
                or self.chains[self.residue_store.i_chain[last_res_idx]] != chain
                or self.residue_store.i_structure[last_res_idx] != self.i_structure
            ):
                is_new_res = True

        if is_new_res:
            self.add_residue(i_atom, res_num, ins_code, chain, res_type)

        i_res = self.get_residue_count() - 1
        self.residue_store.atom_count[i_res] += 1
        self.atom_store.i_res[i_atom] = i_res

        self.atom_store.model[i_atom] = model

    def add_residue(
        self,
        i_first_atom_in_res: int,
        res_num: int,
        ins_code: str,
        chain: str,
        res_type: str,
    ):
        i_res = self.get_residue_count()
        self.residue_store.increment()

        res_id = f"{self.structure_id}:{self.count_res_added}["
        if chain:
            res_id += chain
        res_id += str(res_num) + ins_code.strip()
        res_id += "]"
        self.count_res_added += 1

        self.res_ids.append(res_id)

        i_chain = get_index_from_values(self.chains, chain)
        self.residue_store.i_chain[i_res] = i_chain

        self.residue_store.res_num[i_res] = res_num
        self.residue_store.ins_code[i_res] = char_to_int(ins_code)

        self.residue_store.i_res_type[i_res] = get_index_from_values(
            self.res_types, res_type
        )

        self.residue_store.atom_offset[i_res] = i_first_atom_in_res
        self.residue_store.atom_count[i_res] = 0

        self.residue_store.i_structure[i_res] = self.i_structure

    def find_first_residue(self, chain: str, res_num: int, pdb_id: str = None):
        # This would use actual residue proxy in real implementation
        for i_res in range(self.get_residue_count()):
            if pdb_id:
                if self.structure_ids[self.residue_store.i_structure[i_res]] != pdb_id:
                    continue
            if (
                self.chains[self.residue_store.i_chain[i_res]] == chain
                and self.residue_store.res_num[i_res] == res_num
            ):
                return i_res  # Would return residue proxy in real implementation
        return None

    def find_residue_indices(
        self, i_structure: int, chain: str, res_num: int
    ) -> List[int]:
        result = []
        for i_res in range(self.get_residue_count()):
            if (
                self.residue_store.i_structure[i_res] == i_structure
                and self.chains[self.residue_store.i_chain[i_res]] == chain
                and self.residue_store.res_num[i_res] == res_num
            ):
                result.append(i_res)
        return result

    def get_atom_proxy(self, i_atom: Optional[int] = None):
        return AtomProxy(self, i_atom)

    def get_atom_count(self) -> int:
        return self.atom_store.count

    def get_residue_proxy(self, i_res: Optional[int] = None):
        return ResidueProxy(self, i_res)

    def get_residue_count(self) -> int:
        return self.residue_store.count

    def are_close_residues(self, i_res0: int, i_res1: int) -> bool:
        atom0 = self.get_atom_proxy()
        atom1 = self.get_atom_proxy()
        res0 = self.get_residue_proxy(i_res0)
        pos0 = atom0.load(res0.i_atom).pos.copy()
        atom_indices0 = res0.get_atom_indices()
        res1 = self.get_residue_proxy(i_res1)
        pos1 = atom1.load(res1.i_atom).pos.copy()
        atom_indices1 = res1.get_atom_indices()

        if pos_distance(pos0, pos1) > 17:
            return False

        for i_atom0 in atom_indices0:
            for i_atom1 in atom_indices1:
                if pos_distance(atom0.load(i_atom0).pos, atom1.load(i_atom1).pos) < 4:
                    return True

        return False

    def assign_residue_properties(self, i_structure: int):
        # This would be implemented with actual data and proxy classes
        pass

    def get_i_atom_closest(
        self, pos: Tuple[float, float, float], atom_indices: List[int]
    ) -> int:
        i_atom_closest = None
        min_d = 1e6

        for i_atom in atom_indices:
            if i_atom_closest is None:
                i_atom_closest = i_atom
            else:
                # Calculate distance (simplified - would use actual 3D vector math)
                atom_pos = (
                    self.atom_store.x[i_atom],
                    self.atom_store.y[i_atom],
                    self.atom_store.z[i_atom],
                )
                d = pos_distance(pos, atom_pos)
                if d < min_d:
                    i_atom_closest = i_atom
                    min_d = d

        return i_atom_closest

    def get_i_atom_at_position(self, pos: Tuple[float, float, float]) -> int:
        atom_indices = list(range(self.get_atom_count()))
        i_atom_closest = None
        min_d = 1e6

        for i_atom in atom_indices:
            if i_atom_closest is None:
                i_atom_closest = i_atom
            else:
                atom_pos = (
                    self.atom_store.x[i_atom],
                    self.atom_store.y[i_atom],
                    self.atom_store.z[i_atom],
                )
                d = pos_distance(pos, atom_pos)
                if d < min_d:
                    i_atom_closest = i_atom
                    min_d = d

        if min_d < 0.1:
            return i_atom_closest
        else:
            return -1

    def get_center(self, atom_indices: Optional[List[int]] = None) -> Vector3d:
        if atom_indices is None:
            atom_indices = list(range(self.get_atom_count()))

        sum_x = sum(self.atom_store.x[i] for i in atom_indices)
        sum_y = sum(self.atom_store.y[i] for i in atom_indices)
        sum_z = sum(self.atom_store.z[i] for i in atom_indices)

        n = len(atom_indices)
        return Vector3d(sum_x / n, sum_y / n, sum_z / n)

    def calc_max_length(self, atom_indices: List[int] = None) -> float:
        if atom_indices is None:
            atom_indices = list(range(self.get_atom_count()))

        if not atom_indices:
            return 0.0

        # Calculate min/max for each dimension
        min_coords = [float("inf")] * 3
        max_coords = [float("-inf")] * 3

        for i_atom in atom_indices:
            coords = [
                self.atom_store.x[i_atom],
                self.atom_store.y[i_atom],
                self.atom_store.z[i_atom],
            ]

            for dim in range(3):
                if coords[dim] < min_coords[dim]:
                    min_coords[dim] = coords[dim]
                if coords[dim] > max_coords[dim]:
                    max_coords[dim] = coords[dim]

        spans = [max_coords[i] - min_coords[i] for i in range(3)]
        return max(spans)

    def get_extent_from_center(
        self, center: Vector3d, atom_indices: Optional[List[int]] = None
    ) -> float:
        """Get maximum distance from center to any atom in soup."""
        if atom_indices is None:
            atom_indices = list(range(self.get_atom_count()))

        if not atom_indices:
            return 0.0

        max_dist = 0.0
        atom_proxy = self.get_atom_proxy()

        for i_atom in atom_indices:
            atom_proxy.load(i_atom)
            dist = pos_distance(center, atom_proxy.pos)
            if dist > max_dist:
                max_dist = dist

        return 2.0 * max_dist + 4

    def get_atoms_of_chain_containing_residue(self, i_res: int) -> List[int]:
        residue = self.get_residue_proxy(i_res)
        i_structure = residue.i_structure
        chain = residue.chain
        atom_indices = []

        for i in range(self.get_residue_count()):
            residue.i_res = i
            if residue.i_structure == i_structure and residue.chain == chain:
                atom_indices.append(residue.i_atom)

        return atom_indices

    def get_neighbours(self, i_res: int) -> List[int]:
        indices = [i_res]
        for j_res in range(self.get_residue_count()):
            if self.are_close_residues(j_res, i_res):
                indices.append(j_res)
        return indices

    def is_res_close_to_point(self, i_res0: int, pos1) -> bool:
        atom0 = self.get_atom_proxy()
        res0 = self.get_residue_proxy(i_res0)
        pos0 = atom0.load(res0.i_atom).pos.copy()

        if pos_distance(pos0, pos1) > 17:
            return False

        atom_indices0 = res0.get_atom_indices()
        for i_atom0 in atom_indices0:
            if pos_distance(atom0.load(i_atom0).pos, pos1) < 5:
                return True

        return False

    def get_neighbours_of_point(self, pos) -> List[int]:
        indices = []
        for j_res in range(self.get_residue_count()):
            if self.is_res_close_to_point(j_res, pos):
                indices.append(j_res)
        return indices

    def find_atom_in_soup(self, chain, res_num, atom_type) -> AtomProxy:
        """Find atom in soup by chain, residue number, and atom type"""
        if chain == "":
            chain = " "

        atom_proxy = self.get_atom_proxy()
        res_proxy = self.get_residue_proxy()

        for i_atom in range(self.get_atom_count()):
            atom_proxy.load(i_atom)
            res_proxy.load(atom_proxy.i_res)

            if (
                res_proxy.chain == chain
                and res_proxy.res_num == res_num
                and atom_proxy.atom_type == atom_type
            ):
                return atom_proxy

        raise ValueError(
            "Can't find '%s' atom %s of res %d" % (chain, atom_type, res_num)
        )

    def set_atom_bfactors(self, bfactors: List[float]):
        """Set atom ASA values as B-factors in the soup."""
        atom_proxy = self.get_atom_proxy()
        for i_atom in range(self.get_atom_count()):
            atom_proxy.load(i_atom)
            atom_proxy.bfactor = bfactors[i_atom]

    def get_atom_indices(
        self, atom_indices: Optional[Iterable[int]] = None, skip_waters: bool = False
    ) -> List[int]:
        if atom_indices is None:
            atom_indices = list(range(self.get_atom_count()))

        if skip_waters:
            atom_proxy = self.get_atom_proxy()
            residue_proxy = self.get_residue_proxy()

            def res_type(i_atom):
                return residue_proxy.load(atom_proxy.load(i_atom).i_res).res_type

            atom_indices = [
                i_atom for i_atom in atom_indices if res_type(i_atom) != "HOH"
            ]
            n_water = self.get_atom_count() - len(atom_indices)
            print(f"Skipping {n_water} water atoms -> {len(atom_indices)} atoms")

        return atom_indices

