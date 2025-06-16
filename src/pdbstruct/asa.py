#!/usr/bin/env python

import math
import sys
from typing import List, Tuple

import tqdm

from .parse import load_soup, write_pdb
from .soup import Soup
from .spacehash import SpaceHash, vertex_diff_sq

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
    vertices, radii, probe: float = 1.4, n_sphere_point: int = 960
) -> List[float]:
    """
    Returns list of accessible surface areas of the atoms,
    using the probe and atom radius to define the surface.
    """
    spacehash = SpaceHash(vertices)
    sphere_points = generate_sphere_points(n_sphere_point)
    test_point = [0.0, 0.0, 0.0]

    areas = []
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
                if vertex_diff_sq(vertex_j, test_point) < radius_j * radius_j:
                    i_neighbour_start = i_neighbour
                    is_accessible = False
                    break

            if is_accessible:
                n_accessible_point += 1

        const = 4.0 * math.pi / len(sphere_points)
        area = const * n_accessible_point * radius_i * radius_i
        areas.append(area)

    return areas


def calculate_asa_from_soup(
    soup: Soup, probe: float = 1.4, n_sphere_point: int = 960
) -> List[float]:
    """
    Returns list of accessible surface areas of the atoms,
    using the probe and atom radius to define the surface.
    """
    vertices = []
    radii = []
    atom_proxy = soup.get_atom_proxy()
    for i_atom in range(soup.get_atom_count()):
        atom_proxy.load(i_atom)
        vertices.append(atom_proxy.pos.tuple())
        radii.append(atom_proxy.radius)

    return calculate_asa_from_vertices_and_radii(vertices, radii, probe, n_sphere_point)


def calculate_residue_asas(soup: Soup, probe: float = 1.4) -> List[float]:
    """Calculate accessible surface area for each residue."""
    atom_asas = calculate_asa_from_soup(soup, probe)
    residue_asas = []

    residue_proxy = soup.get_residue_proxy()

    for i_res in range(soup.get_residue_count()):
        residue_proxy.load(i_res)
        atom_indices = residue_proxy.get_atom_indices()
        residue_asa = sum(atom_asas[i_atom] for i_atom in atom_indices)
        residue_asas.append(residue_asa)

    return residue_asas


def calculate_fraction_buried(soup: Soup) -> List[float]:
    """Calculate fraction of each residue that is buried."""
    residue_asas = calculate_residue_asas(soup)
    residue_proxy = soup.get_residue_proxy()

    fractions = []
    for i_res in range(soup.get_residue_count()):
        residue_proxy.load(i_res)
        res_type = residue_proxy.res_type

        if res_type in unfolded_ref_asa:
            unfolded_asa = unfolded_ref_asa[res_type]
            fraction = residue_asas[i_res] / unfolded_asa
        else:
            # If residue type not found, assume fully exposed
            fraction = 1.0

        fractions.append(fraction)

    return fractions


def calc_asa(input_file, n_sphere):
    soup = load_soup(input_file)

    if soup.is_empty():
        print("Error: No atoms found in input file")
        sys.exit(1)

    print("Calculating ASA of atoms")
    atom_asas = calculate_asa_from_soup(soup, probe=1.4, n_sphere_point=n_sphere)
    total_asa = sum(atom_asas)
    print(f"Total ASA: {total_asa:.1f} Å²")

    base_name = input_file.replace(".pdb", "").replace(".cif", "")
    output_file = f"{base_name}-asa.pdb"
    soup.set_atom_bfactors(atom_asas)
    write_pdb(soup, output_file)
    print(f"Wrote ASA as atom bfactors in {output_file}")


usage = """
Copyright (c) 2007 Bosco Ho

Calculates the total Accessible Surface Area (ASA) of atoms in a 
PDB file. 

Usage: hollow-asa.py [-n n_sphere] in_pdb [out_pdb]

- out_pdb    PDB file in which the atomic ASA values are written 
           to the b-factor column.

-n n_sphere  number of points used in generating the spherical
           dot-density for the calculation (default=960). The 
           more points, the more accurate (but slower) the 
           calculation. 
"""


def main():
    import getopt
    import sys

    try:
        opts, args = getopt.getopt(sys.argv[1:], "n:")
    except getopt.GetoptError as e:
        print(f"Error: {e}")
        print(usage)
        sys.exit(1)

    if len(args) == 0:
        print(usage)
        sys.exit(1)

    # Parse command line options
    n_sphere = 960
    for o, a in opts:
        if "-n" in o:
            try:
                n_sphere = int(a)
                print(f"Points on sphere: {n_sphere}")
            except ValueError:
                print("Error: -n option requires an integer value")
                sys.exit(1)

    input_file = args[0]

    try:
        calc_asa(input_file, n_sphere)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
