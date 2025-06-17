import os
import re
from typing import Iterable, List, Optional

from .soup import Soup


def delete_numbers(text: str) -> str:
    """Remove the first occurrence of digits from text."""
    return re.sub(r"\d+", "", text, count=1)


def remove_quotes(s: str) -> str:
    """Remove surrounding quotes from a string."""
    n = len(s)
    if n >= 2:
        if (s[0] == '"' and s[n - 1] == '"') or (s[0] == "'" and s[n - 1] == "'"):
            return s[1 : n - 1]
    return s


def pad_atom_type(in_atom_type):
    atom_type = in_atom_type
    if len(atom_type) == 1:
        atom_type = " %s  " % atom_type
    elif len(atom_type) == 2:
        atom_type = " %s " % atom_type
    elif len(atom_type) == 3:
        if atom_type[0].isdigit():
            atom_type = "%s " % atom_type
        else:
            atom_type = " %s" % atom_type
    return atom_type


def strip_left_if_too_long(s, max_length):
    """
    Returns the string, stripped from the left if it exceeds max_length.
    Keeps only the rightmost max_length characters.
    """
    if len(s) > max_length:
        return s[-max_length:]
    return s


def strip_lines(pdb_txt, tag_func):
    new_lines = []
    for line in pdb_txt.splitlines():
        if tag_func(line):
            continue
        new_lines.append(line)
    return "\n".join(new_lines)


def strip_hydrogens(pdb_txt):
    new_lines = []
    for line in pdb_txt.splitlines():
        if line.startswith("ATOM"):
            atom_type = line[12:16]
            if atom_type[0] == "H" or atom_type[1] == "H":
                continue
        new_lines.append(line)
    return "\n".join(new_lines)


def strip_other_nmr_models(pdb_txt):
    new_lines = []
    for line in pdb_txt.splitlines():
        new_lines.append(line)
        if line.startswith("ENDMDL"):
            break
    return "\n".join(new_lines)


def strip_alternative_atoms(pdb_txt):
    new_lines = []
    for line in pdb_txt.splitlines():
        new_line = line
        if line.startswith("ATOM"):
            alt_loc = line[16]
            if alt_loc not in [" "]:
                if alt_loc in ["A", "a"]:
                    new_line = line[:16] + " " + line[17:]
                else:
                    continue
        new_lines.append(new_line)
    return "\n".join(new_lines)


class PdbParser:
    def __init__(self, soup: Soup, scrub: bool = False, skip_water: bool=False) -> None:
        self.soup: Soup = soup
        self.scrub = scrub
        self.skip_water = skip_water
        self.has_secondary_structure = False
        self.error = ""

    def is_atom_line(self, line: str) -> bool:
        return line.startswith("ATOM") or line.startswith("HETATM")

    def is_nmr(self, lines: List[str]) -> bool:
        for line in lines:
            if line.startswith("EXPDTA"):
                if "NMR" in line:
                    return True
        return False

    def parse_atom_lines(self, pdb_lines: List[str], model: int) -> None:
        self.soup.count_res_added = 0

        for i_line, line in enumerate(pdb_lines):
            if self.is_atom_line(line):
                try:
                    atom_type = line[12:16].strip()
                    alt = line[16:17].strip()
                    res_type = line[17:20].strip()
                    chain = line[21]
                    res_num = int(line[22:26])
                    ins_code = line[26:27]
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    occupancy = float(line[54:60])
                    bfactor = float(line[60:66])
                    elem = line[76:78].strip()
                except (ValueError, IndexError):
                    self.error = f"line {i_line}"
                    print(f'parse_atom_lines: "{line}"')
                    continue

                if elem == "":
                    elem = delete_numbers(atom_type.strip())[:1]

                if self.scrub:
                    if alt not in ["", " ", "A", "a"]:
                        continue
                    if model > 1:
                        continue

                if self.skip_water:
                    if res_type == "HOH":
                        continue

                self.soup.add_atom(
                    x=x,
                    y=y,
                    z=z,
                    bfactor=bfactor,
                    alt=alt,
                    atom_type=atom_type,
                    elem=elem,
                    res_type=res_type,
                    res_num=res_num,
                    ins_code=ins_code,
                    chain=chain,
                    occupancy=occupancy,
                    model=model,
                )

    def parse_secondary_structure_lines(self, pdb_lines: List[str]) -> None:
        self.soup.assign_residue_properties(self.soup.i_structure)
        residue = self.soup.get_residue_proxy()

        for i_line, line in enumerate(pdb_lines):
            if line.startswith("HELIX"):
                self.has_secondary_structure = True
                chain = line[19:20]
                res_num_start = int(line[21:25])
                res_num_end = int(line[33:37])

                for i_res in self.soup.find_residue_indices(
                    self.soup.i_structure, chain, res_num_start
                ):
                    residue.i_res = i_res
                    while (
                        residue.i_res < self.soup.get_residue_count()
                        and residue.res_num <= res_num_end
                        and chain == residue.chain
                    ):
                        residue.ss = "H"
                        residue.i_res = residue.i_res + 1

            elif line.startswith("SHEET"):
                self.has_secondary_structure = True
                chain = line[21:22]
                res_num_start = int(line[22:26])
                res_num_end = int(line[33:37])

                for i_res in self.soup.find_residue_indices(
                    self.soup.i_structure, chain, res_num_start
                ):
                    residue.i_res = i_res
                    while (
                        residue.i_res < self.soup.get_residue_count()
                        and residue.res_num <= res_num_end
                        and chain == residue.chain
                    ):
                        residue.ss = "E"
                        residue.i_res = residue.i_res + 1

    def parse_title(self, lines: List[str]) -> str:
        result = ""
        for line in lines:
            if line[:5] == "TITLE":
                result += line[10:]
        return result

    def parse_text(self, text: str, pdb_id: str) -> None:
        lines = text.split("\n")
        # Handle both \n and \r\n line endings
        lines = [line.rstrip("\r") for line in lines]

        if len(lines) == 0:
            self.parsing_error = "No atom lines"
            return

        title = self.parse_title(lines)
        is_nmr = self.is_nmr(lines)

        lines_list = [[]]
        i_model = 0

        for line in lines:
            if self.is_atom_line(line):
                lines_list[i_model].append(line)
            elif line.startswith("END"):
                if is_nmr:
                    break
                lines_list.append([])
                i_model += 1

        # Remove empty models at the end
        while i_model >= 0 and len(lines_list[i_model]) == 0:
            lines_list.pop()
            i_model -= 1

        n_model = len(lines_list)
        for i_model in range(n_model):
            structure_id = pdb_id
            if n_model > 1:
                structure_id = f"{structure_id}[{i_model + 1}]"

            base_structure_id = structure_id
            i_clash = 1
            while structure_id in self.soup.structure_ids:
                structure_id = f"{base_structure_id}[{i_clash}]"
                i_clash += 1

            self.soup.push_structure_id(structure_id, title)
            model_num = i_model + 1
            lines = lines_list[i_model]
            self.parse_atom_lines(lines, model_num)
            self.parse_secondary_structure_lines(lines)


class CifParser:
    def __init__(self, soup: Soup, scrub: bool = False, skip_water: bool = False):
        self.soup: Soup = soup
        self.scrub = scrub
        self.skip_water = skip_water
        self.has_secondary_structure = False
        self.error = ""

    def is_atom_line(self, line: str) -> bool:
        return line.startswith("ATOM") or line.startswith("HETATM")

    def parse_atom_lines(self, pdb_lines: List[str]) -> None:
        """Parse atom lines from CIF format."""
        next_res_num = None
        last_chain = None
        last_entity = None

        for i_line, line in enumerate(pdb_lines):
            if self.is_atom_line(line):
                tokens = re.split(r"[ ,]+", line)
                try:
                    elem = tokens[2]
                    atom_type = remove_quotes(tokens[3])
                    alt = "" if tokens[4] == "." else tokens[4]
                    res_type = tokens[5]
                    chain = tokens[6]
                    entity = tokens[7]
                    ins_code = " " if tokens[9] == "?" else tokens[9]
                    x = float(tokens[10])
                    y = float(tokens[11])
                    z = float(tokens[12])
                    occupancy = float(tokens[13])
                    bfactor = float(tokens[14])

                    # res_num
                    if tokens[8] == ".":
                        is_same_chain_and_entity = (
                            chain == last_chain and entity == last_entity
                        )
                        if not is_same_chain_and_entity or res_type == "HOH":
                            if next_res_num is None:
                                next_res_num = 1
                            else:
                                next_res_num += 1
                            last_chain = chain
                            last_entity = entity
                        res_num = next_res_num
                    else:
                        res_num = int(tokens[16])
                        last_chain = chain
                        last_entity = entity
                        next_res_num = res_num + 1

                    model = int(tokens[20])

                except (ValueError, IndexError) as e:
                    self.error = f"line {i_line}"
                    print(f'parse_atom_lines {e}: "{line}"')
                    continue

                if elem == "":
                    elem = delete_numbers(atom_type.strip())[:1]

                if self.scrub:
                    if alt not in ["", " ", "A", "a"]:
                        continue
                    if model > 1:
                        continue

                if self.skip_water:
                    if res_type == "HOH":
                        continue

                self.soup.add_atom(
                    x=x,
                    y=y,
                    z=z,
                    bfactor=bfactor,
                    alt=alt,
                    atom_type=atom_type,
                    elem=elem,
                    res_type=res_type,
                    res_num=res_num,
                    ins_code=ins_code,
                    chain=chain,
                    occupancy=occupancy,
                    model=model,
                )

    def parse_secondary_structure_lines(self, pdb_lines: List[str]) -> None:
        self.has_secondary_structure = False
        self.soup.assign_residue_properties(self.soup.i_structure)
        self.parse_helix_lines(pdb_lines)
        self.parse_sheet_lines(pdb_lines)

    def parse_helix_lines(self, pdb_lines: List[str]) -> None:
        residue = self.soup.get_residue_proxy()
        is_helix_loop = False

        for i_line, line in enumerate(pdb_lines):
            if not is_helix_loop:
                if line.startswith("_struct_conf.pdbx_PDB_helix_id"):
                    is_helix_loop = True
                continue

            if line.startswith("#"):
                break

            if not line.startswith("_struct_conf"):
                self.has_secondary_structure = True
                tokens = re.split(r"[ ,]+", line)
                chain = tokens[4]
                res_num_start = int(tokens[5])
                res_num_end = int(tokens[9])

                for i_res in self.soup.find_residue_indices(
                    self.soup.i_structure, chain, res_num_start
                ):
                    residue.i_res = i_res
                    while residue.res_num <= res_num_end and chain == residue.chain:
                        residue.ss = "H"
                        residue.i_res = residue.i_res + 1

    def parse_sheet_lines(self, pdb_lines: List[str]) -> None:
        residue = self.soup.get_residue_proxy()
        is_sheet_loop = False

        for i_line, line in enumerate(pdb_lines):
            if not is_sheet_loop:
                if line.startswith("_struct_sheet_range.sheet_id"):
                    is_sheet_loop = True
                continue

            if line.startswith("#"):
                break

            if not line.startswith("_struct"):
                self.has_secondary_structure = True
                tokens = re.split(r"[ ,]+", line)
                chain = tokens[3]
                res_num_start = int(tokens[4])
                res_num_end = int(tokens[8])

                for i_res in self.soup.find_residue_indices(
                    self.soup.i_structure, chain, res_num_start
                ):
                    residue.i_res = i_res
                    while residue.res_num <= res_num_end and chain == residue.chain:
                        residue.ss = "E"
                        residue.i_res = residue.i_res + 1

    def parse_title(self, lines: List[str]) -> str:
        for i, line in enumerate(lines):
            if line.startswith("_struct.title"):
                rest = line.replace("_struct.title", "").strip()
                if rest:
                    return remove_quotes(rest)
            if i > 0:
                prev_line = lines[i - 1]
                if prev_line.startswith("_struct.title"):
                    return remove_quotes(line.strip())
        return ""

    def parse_text(self, text: str, pdb_id: str) -> None:
        lines = text.split("\n")
        lines = [line.rstrip("\r") for line in lines]

        if len(lines) == 0:
            self.parsing_error = "No atom lines"
            return

        self.soup.push_structure_id(pdb_id, self.parse_title(lines))
        self.parse_atom_lines(lines)
        self.parse_secondary_structure_lines(lines)


def load_soup(filename: str, scrub=False, skip_water=False) -> Soup:
    """Load structure from PDB or CIF file."""
    soup = Soup()

    # Determine file type and parse accordingly
    ext = os.path.splitext(filename)[1].lower()
    if ext in [".cif", ".mmcif"]:
        parser = CifParser(soup, scrub, skip_water)
    elif ext in [".pdb", ".pdbx"]:
        parser = PdbParser(soup, scrub, skip_water)
    else:
        raise ValueError(f"Unknown file type: {ext}")

    # Read file content
    try:
        with open(filename, "r") as f:
            content = f.read()
    except IOError as e:
        raise IOError(f"Could not read file {filename}: {e}")

    pdb_id = os.path.splitext(os.path.basename(filename))[0]

    parser.parse_text(content, pdb_id)

    if parser.error:
        print(f"Warning: Parser error - {parser.error}")

    return soup


def add_suffix_to_basename(filename: str, suffix: str) -> str:
    """
    Add a suffix to the basename of a filename while preserving directory and extension.

    Args:
        filename: Input filename (can include path)
        suffix: Suffix to add to basename (e.g., "-extra")

    Returns:
        Modified filename with suffix added to basename

    Examples:
        add_suffix_to_basename("protein.pdb", "-extra") -> "protein-extra.pdb"
        add_suffix_to_basename("/path/to/protein.pdb", "-extra") -> "/path/to/protein-extra.pdb"
        add_suffix_to_basename("protein", "-extra") -> "protein-extra"
    """
    import os

    # Split into directory, basename, and extension
    dir_path = os.path.dirname(filename)
    base_name = os.path.basename(filename)

    # Split basename into name and extension
    name, ext = os.path.splitext(base_name)

    # Add suffix to name
    new_name = name + suffix

    # Reconstruct the full path
    new_basename = new_name + ext

    if dir_path:
        return os.path.join(dir_path, new_basename)
    else:
        return new_basename


def write_pdb(soup: Soup, filename: str, atom_indices: Optional[Iterable[int]] = None) -> None:
    atom_proxy = soup.get_atom_proxy()
    residue_proxy = soup.get_residue_proxy()

    with open(filename, "w") as f:
        if soup.title:
            f.write(f"TITLE     {soup.title}\n")

        for i_atom in soup.get_atom_indices(atom_indices):
            atom_proxy.load(i_atom)
            residue_proxy.load(atom_proxy.i_res)

            atom_counter = strip_left_if_too_long(str(i_atom + 1), 5)
            res_num = strip_left_if_too_long(str(residue_proxy.res_num), 4)

            # Format ATOM record according to PDB specification
            # Columns: 1-6: Record name, 7-11: Atom serial, 13-16: Atom name,
            # 17: Alt loc, 18-20: Residue name, 22: Chain, 23-26: Residue seq,
            # 27: Insertion code, 31-38: X, 39-46: Y, 47-54: Z,
            # 55-60: Occupancy, 61-66: B-factor, 77-78: Element
            line = (
                f"{'ATOM':<6}{atom_counter:>5} {pad_atom_type(atom_proxy.atom_type):<4}"
                f"{atom_proxy.alt if atom_proxy.alt else ' ':1}{residue_proxy.res_type:<3} "
                f"{residue_proxy.chain:1}{res_num:>4}"
                f"{residue_proxy.ins_code if residue_proxy.ins_code else ' ':1}   "
                f"{atom_proxy.pos.x:8.3f}{atom_proxy.pos.y:8.3f}{atom_proxy.pos.z:8.3f}"
                f"{atom_proxy.occupancy:6.2f}{atom_proxy.bfactor:6.2f}          {atom_proxy.elem:>2}\n"
            )
            f.write(line)

        f.write("END\n")


def write_cif(soup: Soup, filename: str, atom_indices: Optional[Iterable[int]] = None):
    """Write structure data to CIF format file."""
    atom_proxy = soup.get_atom_proxy()
    residue_proxy = soup.get_residue_proxy()

    with open(filename, "w") as f:
        # Write header information
        f.write("data_structure\n")
        f.write("#\n")

        if soup.title:
            f.write(f"_struct.title '{soup.title}'\n")
        f.write("#\n")

        # Write atom site loop header
        f.write("loop_\n")
        f.write("_atom_site.group_PDB\n")
        f.write("_atom_site.id\n")
        f.write("_atom_site.type_symbol\n")
        f.write("_atom_site.label_atom_id\n")
        f.write("_atom_site.label_alt_id\n")
        f.write("_atom_site.label_comp_id\n")
        f.write("_atom_site.label_asym_id\n")
        f.write("_atom_site.label_entity_id\n")
        f.write("_atom_site.label_seq_id\n")
        f.write("_atom_site.pdbx_PDB_ins_code\n")
        f.write("_atom_site.Cartn_x\n")
        f.write("_atom_site.Cartn_y\n")
        f.write("_atom_site.Cartn_z\n")
        f.write("_atom_site.occupancy\n")
        f.write("_atom_site.B_iso_or_equiv\n")
        f.write("_atom_site.pdbx_formal_charge\n")
        f.write("_atom_site.auth_seq_id\n")
        f.write("_atom_site.auth_comp_id\n")
        f.write("_atom_site.auth_asym_id\n")
        f.write("_atom_site.auth_atom_id\n")
        f.write("_atom_site.pdbx_PDB_model_num\n")

        # Write atom data
        entity_id = 1
        current_chain = None

        for i_atom in soup.get_atom_indices(atom_indices):
            atom_proxy.load(i_atom)
            residue_proxy.load(atom_proxy.i_res)

            # Determine if this is ATOM or HETATM
            group_pdb = (
                "ATOM"
                if residue_proxy.res_type
                in [
                    "ALA",
                    "ARG",
                    "ASN",
                    "ASP",
                    "CYS",
                    "GLN",
                    "GLU",
                    "GLY",
                    "HIS",
                    "ILE",
                    "LEU",
                    "LYS",
                    "MET",
                    "PHE",
                    "PRO",
                    "SER",
                    "THR",
                    "TRP",
                    "TYR",
                    "VAL",
                ]
                else "HETATM"
            )

            # Update entity_id when chain changes
            if current_chain != residue_proxy.chain:
                if current_chain is not None:
                    entity_id += 1
                current_chain = residue_proxy.chain

            # Format atom record
            atom_id = i_atom + 1
            alt_id = atom_proxy.alt if atom_proxy.alt else "."
            ins_code = residue_proxy.ins_code if residue_proxy.ins_code else "?"

            # Write the atom line
            f.write(f"{group_pdb:<6} ")
            f.write(f"{atom_id:<6} ")
            f.write(f"{atom_proxy.elem:<2} ")
            f.write(f"{atom_proxy.atom_type:<4} ")
            f.write(f"{alt_id:<1} ")
            f.write(f"{residue_proxy.res_type:<3} ")
            f.write(f"{residue_proxy.chain:<1} ")
            f.write(f"{entity_id:<1} ")
            f.write(f"{residue_proxy.res_num:<4} ")
            f.write(f"{ins_code:<1} ")
            f.write(f"{atom_proxy.pos.x:8.3f} ")
            f.write(f"{atom_proxy.pos.y:8.3f} ")
            f.write(f"{atom_proxy.pos.z:8.3f} ")
            f.write(f"{atom_proxy.occupancy:6.2f} ")
            f.write(f"{atom_proxy.bfactor:6.2f} ")
            f.write("? ")  # pdbx_formal_charge
            f.write(f"{residue_proxy.res_num:<4} ")
            f.write(f"{residue_proxy.res_type:<3} ")
            f.write(f"{residue_proxy.chain:<1} ")
            f.write(f"{atom_proxy.atom_type:<4} ")
            f.write("1")  # pdbx_PDB_model_num
            f.write("\n")

        f.write("#\n")


def write_soup(soup: Soup, filename: str, atom_indices: Optional[Iterable[int]] = None):
    """
    Write structure data to file, automatically choosing format based on extension.

    Args:
        soup: Soup object containing structure data
        filename: Output filename (extension determines format)

    Raises:
        ValueError: If file extension is not supported
    """
    ext = os.path.splitext(filename)[1].lower()

    if ext == ".cif":
        write_cif(soup, filename, atom_indices)
    elif ext in [".pdb", ".ent"]:
        write_pdb(soup, filename, atom_indices)
    else:
        raise ValueError(f"Unsupported file extension '{ext}'. Use .pdb, .ent, or .cif")
