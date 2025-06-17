#!/usr/bin/env python

import math
import sys
from typing import List, Optional, Set, Tuple

import tqdm

from .parse import add_suffix_to_basename, load_soup, write_soup
from .soup import Soup
from .spacehash import SpaceHash
from .vector3d import pos_distance_sq

__doc__ = """
Routines to calculate the Accessible Surface Area of a set of atoms.
The algorithm is adapted from the Rose lab's chasa.py, which uses
the dot density technique found in:

Shrake, A., and J. A. Rupley. "Environment and Exposure to Solvent
of Protein Atoms. Lysozyme and Insulin." JMB (1973) 79:351-371.
"""


# Reference ASA values for unfolded proteins
# taken from fragments of length 17 from
# Creamer et al. (1995) Biochemistry: 36:2832
unfolded_ref_asa = {
    "ALA": 19.8 + 46.6,
    "ARG": 17.1 + 156.9,
    "ASN": 17.6 + 84.5,
    "ASP": 18.1 + 79.2,
    "CYS": 18.2 + 62.9,
    "GLN": 17.2 + 105.0,
    "GLU": 17.9 + 102.8,
    "GLY": 54.6 + 0.0,
    "HIS": 14.9 + 103.9,
    "ILE": 15.2 + 100.1,
    "LEU": 14.7 + 101.4,
    "LYS": 18.3 + 142.5,
    "MET": 16.7 + 105.3,
    "PHE": 15.3 + 118.7,
    "PRO": 18.9 + 83.5,
    "SER": 23.8 + 59.7,
    "THR": 18.6 + 77.3,
    "TRP": 15.1 + 154.7,
    "TYR": 17.7 + 131.0,
    "VAL": 15.9 + 81.8,
}


def generate_sphere_points(n: int) -> List[Tuple[float, float, float]]:
    """
    Returns list of 3d coordinates of points on a sphere using the
    Golden Section Spiral algorithm.
    """
    if n == 0:
        return []
    points = []
    inc = math.pi * (3 - math.sqrt(5))
    offset = 2 / float(n)
    for k in range(int(n)):
        y = k * offset - 1 + (offset / 2.0)
        r = math.sqrt(1 - y * y)
        phi = k * inc
        points.append((math.cos(phi) * r, y, math.sin(phi) * r))
    return points


def reorder_range(n, i_start):
    return list(range(i_start, n)) + list(range(i_start))


def calculate_asa_from_vertices_and_radii(
    vertices,
    radii,
    probe: float = 1.4,
    n_sphere_point: int = 960,
) -> List[float]:
    """
    Returns list of accessible surface areas of the atoms,
    using the probe and atom radius to define the surface.

    Args:
        vertices: List of vertex coordinates
        radii: List of atomic radii
        probe: Probe radius for surface calculation
        n_sphere_point: Number of sphere points for calculation
    """
    spacehash = SpaceHash(vertices)
    sphere_points = generate_sphere_points(n_sphere_point)
    test_point = [0.0, 0.0, 0.0]

    # Initialize areas list with zeros for all vertices
    areas = [0.0] * len(vertices)

    for i_vertex in tqdm.trange(len(vertices)):
        neighbor_vertex_indices = spacehash.find_connected_vertex_indices(
            radii, probe, i_vertex
        )
        n_neighbor = len(neighbor_vertex_indices)
        i_neighbour_start = 0
        radius_i = probe + radii[i_vertex]
        vertex_i = spacehash.vertices[i_vertex]
        n_accessible_point = 0

        for sphere_point in sphere_points:
            is_accessible = True

            test_point[0] = sphere_point[0] * radius_i + vertex_i[0]
            test_point[1] = sphere_point[1] * radius_i + vertex_i[1]
            test_point[2] = sphere_point[2] * radius_i + vertex_i[2]

            for i_neighbour in reorder_range(n_neighbor, i_neighbour_start):
                j_vertex = neighbor_vertex_indices[i_neighbour]
                radius_j = radii[j_vertex] + probe
                vertex_j = spacehash.vertices[j_vertex]
                if pos_distance_sq(vertex_j, test_point) < radius_j * radius_j:
                    i_neighbour_start = i_neighbour
                    is_accessible = False
                    break

            if is_accessible:
                n_accessible_point += 1

        fraction = n_accessible_point / len(sphere_points)
        areas[i_vertex] = fraction * 4.0 * math.pi * radius_i * radius_i

    return areas


def calculate_asa_from_soup(
    soup: Soup,
    probe: float = 1.4,
    n_sphere_point: int = 960,
    atom_indices: Optional[List[int]] = None,
) -> List[float]:
    """
    Returns list of accessible surface areas of the atoms,
    using the probe and atom radius to define the surface.

    Args:
        soup: Soup object containing atomic data
        probe: Probe radius for surface calculation
        n_sphere_point: Number of sphere points for calculation
        atom_indices: Set of atom indices to calculate ASA for.
                              If None, calculates for all atoms.
    """
    if atom_indices is None:
        atom_indices = range(soup.get_atom_count())

    vertices = []
    radii = []
    atom_proxy = soup.get_atom_proxy()
    for i_atom in atom_indices:
        atom_proxy.load(i_atom)
        vertices.append(atom_proxy.pos.tuple())
        radii.append(atom_proxy.radius)

    return calculate_asa_from_vertices_and_radii(
        vertices, radii, probe, n_sphere_point
    )


def calculate_residue_asas(
    soup: Soup, probe: float = 1.4, selected_residue_indices: Optional[Set[int]] = None
) -> List[float]:
    """
    Calculate accessible surface area for each residue.

    Args:
        soup: Soup object containing atomic data
        probe: Probe radius for surface calculation
        selected_residue_indices: Set of residue indices to calculate ASA for.
                                 If None, calculates for all residues.
    """
    # Get atom indices for selected residues
    atom_indices = None
    if selected_residue_indices is not None:
        atom_indices = set()
        residue_proxy = soup.get_residue_proxy()
        for i_res in selected_residue_indices:
            if 0 <= i_res < soup.get_residue_count():
                residue_proxy.load(i_res)
                atom_indices = residue_proxy.get_atom_indices()
                atom_indices.update(atom_indices)

    selected_asas = calculate_asa_from_soup(soup, probe, atom_indices=atom_indices)

    if atom_indices is None:
        atom_asas = selected_asas
    else:
        atom_asas = [0.0] * soup.get_atom_count()
        for i_atom, asa in zip(atom_indices, selected_asas):
            atom_asas[i_atom] = asa

    # Initialize residue ASAs with zeros
    residue_asas = [0.0] * soup.get_residue_count()
    residue_proxy = soup.get_residue_proxy()

    # Determine which residues to process
    if selected_residue_indices is None:
        residue_indices_to_process = range(soup.get_residue_count())
    else:
        residue_indices_to_process = [
            i for i in selected_residue_indices if 0 <= i < soup.get_residue_count()
        ]

    for i_res in residue_indices_to_process:
        residue_proxy.load(i_res)
        atom_indices = residue_proxy.get_atom_indices()
        residue_asa = sum(atom_asas[i_atom] for i_atom in atom_indices)
        residue_asas[i_res] = residue_asa

    return residue_asas


def calculate_fraction_buried(
    soup: Soup, selected_residue_indices: Optional[Set[int]] = None
) -> List[float]:
    """
    Calculate fraction of each residue that is buried.

    Args:
        soup: Soup object containing atomic data
        selected_residue_indices: Set of residue indices to calculate for.
                                 If None, calculates for all residues.
    """
    residue_asas = calculate_residue_asas(
        soup, selected_residue_indices=selected_residue_indices
    )
    residue_proxy = soup.get_residue_proxy()

    # Initialize fractions with zeros
    fractions = [0.0] * soup.get_residue_count()

    # Determine which residues to process
    if selected_residue_indices is None:
        residue_indices_to_process = range(soup.get_residue_count())
    else:
        residue_indices_to_process = [
            i for i in selected_residue_indices if 0 <= i < soup.get_residue_count()
        ]

    for i_res in residue_indices_to_process:
        residue_proxy.load(i_res)
        res_type = residue_proxy.res_type

        if res_type in unfolded_ref_asa:
            unfolded_asa = unfolded_ref_asa[res_type]
            fraction = residue_asas[i_res] / unfolded_asa if unfolded_asa > 0 else 0.0
        else:
            # If residue type not found, assume fully exposed
            fraction = 1.0

        fractions[i_res] = fraction

    return fractions


def calc_asa(input_file, n_sphere, skip_waters: bool=False):
    """
    Calculate ASA for atoms in a PDB file.

    Args:
        input_file: Path to input PDB file
        n_sphere: Number of sphere points for calculation
        selected_atom_indices: Set of atom indices to calculate ASA for.
                              If None, calculates for all atoms.
    """
    soup = load_soup(input_file, scrub=True)
    print(
        f"Loaded {soup.get_atom_count()} atoms in {soup.get_residue_count()} residues from `{input_file}`"
    )

    if soup.is_empty():
        print("Error: No atoms found in input file")
        sys.exit(1)

    atom_indices = soup.get_atom_indices(skip_waters=skip_waters)

    print("Calculating ASA of atoms")
    selected_atom_asas = calculate_asa_from_soup(
        soup,
        probe=1.4,
        n_sphere_point=n_sphere,
        atom_indices=atom_indices,
    )

    # Calculate total ASA only for selected atoms or all atoms
    total_asa = sum(selected_atom_asas)

    print(f"Total ASA: {total_asa:.1f} Å²")

    all_atom_asas = [0.0] * soup.get_atom_count()
    for i_atom, asa in zip(atom_indices, selected_atom_asas):
        all_atom_asas[i_atom] = asa
    soup.set_atom_bfactors(all_atom_asas)
    output_file = add_suffix_to_basename(input_file, "-asa")
    write_soup(soup, output_file, atom_indices=atom_indices)
    print(f"Wrote ASA as atom bfactors in {output_file}")
