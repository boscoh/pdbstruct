#!/usr/bin/env python

import math
import os
import sys

import click
import tqdm

from . import asa, vector3d
from .bgrid import BoolGrid
from .util import read_parameters, this_dir
from .parse import add_suffix_to_basename, load_soup, write_soup
from .soup import Soup
from .spacehash import SpaceHash
from .vector3d import pos_distance_sq


defaults = read_parameters(os.path.join(this_dir, "hollow.defaults.txt"))


class HollowGrid:
    def __init__(self, grid_spacing, width, center):
        self.center = center
        self.width = width
        self.spacing = grid_spacing
        self.inv_spacing = 1.0 / self.spacing

        self.n = int(math.ceil(self.width * self.inv_spacing))
        self.half_n = self.n // 2

        self.excluded_grid = BoolGrid(self.n)
        self.is_excluded = self.excluded_grid.is_set
        self.set_excluded = self.excluded_grid.set

        self.drilled_grid = BoolGrid(self.n)
        self.is_drilled = self.drilled_grid.is_set
        self.set_drilled = self.drilled_grid.set

        self.x = [
            self.center.x + (i - self.half_n) * self.spacing for i in range(self.n)
        ]
        self.y = [
            self.center.y + (i - self.half_n) * self.spacing for i in range(self.n)
        ]
        self.z = [
            self.center.z + (i - self.half_n) * self.spacing for i in range(self.n)
        ]

    def is_excluded_or_drilled(self, i, j, k):
        return self.is_excluded(i, j, k) or self.is_drilled(i, j, k)

    def indices(self, pos):
        return (
            (pos.x - self.center.x) * self.inv_spacing + self.half_n,
            (pos.y - self.center.y) * self.inv_spacing + self.half_n,
            (pos.z - self.center.z) * self.inv_spacing + self.half_n,
        )

    def pos(self, i, j, k):
        return vector3d.Vector3d(self.x[i], self.y[j], self.z[k])

    def is_grid_point_near_sphere(self, i, j, k, vertex, r_sq):
        d_x = self.x[i] - vertex[0]
        d_y = self.y[j] - vertex[1]
        d_z = self.z[k] - vertex[2]
        return d_x * d_x + d_y * d_y + d_z * d_z < r_sq

    def int_range(self, low_f, high_f):
        low = max(0, int(math.floor(low_f)))
        high = min(self.n, int(math.ceil(high_f)) + 1)
        return list(range(low, high))

    def exclude_sphere(self, vertex, r):
        r_sq = r * r
        low = vector3d.Vector3d(vertex[0] - r, vertex[1] - r, vertex[2] - r)
        high = vector3d.Vector3d(vertex[0] + r, vertex[1] + r, vertex[2] + r)
        low_i, low_j, low_k = self.indices(low)
        high_i, high_j, high_k = self.indices(high)
        for i in self.int_range(low_i, high_i):
            for j in self.int_range(low_j, high_j):
                for k in self.int_range(low_k, high_k):
                    if not self.is_excluded(i, j, k):
                        if self.is_grid_point_near_sphere(i, j, k, vertex, r_sq):
                            self.set_excluded(i, j, k, True)

    def permutation(self, i, j, k, dim):
        if dim == 0:
            return i, j, k
        if dim == 1:
            return j, k, i
        if dim == 2:
            return k, i, j

    def drill_in_dim(self, is_reversed, i, j, dim):
        drill_range = list(range(self.n))
        if is_reversed:
            drill_range.reverse()
        for k in drill_range:
            a, b, c = self.permutation(i, j, k, dim)
            if self.is_excluded(a, b, c):
                return
            self.set_drilled(a, b, c, True)

    def exclude_edge_to_interior(self):
        for i in range(self.n):
            for j in range(self.n):
                self.drill_in_dim(True, i, j, 0)
                self.drill_in_dim(False, i, j, 0)
                self.drill_in_dim(True, i, j, 1)
                self.drill_in_dim(False, i, j, 1)
                self.drill_in_dim(True, i, j, 2)
                self.drill_in_dim(False, i, j, 2)

    def is_surrounded(self, i, j, k):
        indices_list = [
            (i, j, k),
            (i + 1, j, k),
            (i - 1, j, k),
            (i, j + 1, k),
            (i, j - 1, k),
            (i, j, k - 1),
            (i, j, k + 1),
        ]
        for a, b, c in indices_list:
            if 0 <= a < self.n and 0 <= b < self.n and 0 <= c < self.n:
                if self.is_excluded_or_drilled(a, b, c):
                    return False
        return True

    def exclude_surrounded(self, skip):
        surrounded_grid_points = []
        for i in range(self.n):
            for j in range(self.n):
                for k in range(self.n):
                    if self.is_surrounded(i, j, k):
                        surrounded_grid_points.append([i, j, k])
        for i, j, k in surrounded_grid_points:
            if skip > 0:
                if i % skip == 0 and j % skip == 0 and k % skip == 0:
                    continue
            self.set_excluded(i, j, k, True)

    def exclude_points_in_constraint(self, constraint_fn):
        for i in range(self.n):
            for j in range(self.n):
                for k in range(self.n):
                    if not self.is_excluded_or_drilled(i, j, k):
                        if not constraint_fn(self.pos(i, j, k)):
                            self.set_excluded(i, j, k, True)

    def exclude_vertices(self, vertices, radii, probe):
        for i in tqdm.trange(len(vertices)):
            self.exclude_sphere(vertices[i], radii[i] + probe)

    def exclude_surface(self, vertices, radii, vertex_indices, probe):
        """
        Exclude grid points that are accessible from the molecular surface.

        This method implements a molecular surface accessibility algorithm by rolling
        a probe sphere over the surface of selected atoms. Grid points that fall within
        the probe radius of accessible surface points are marked as excluded.

        The algorithm works by:
        1. For each surface atom, generate test points on a sphere around it
        2. Check if each test point is accessible (not blocked by neighboring atoms)
        3. For accessible test points, exclude all grid points within probe radius

        This effectively creates a "shell" around the accessible molecular surface,
        which is useful for identifying surface cavities and channels while excluding
        bulk solvent regions.

        Args:
            vertices (list): List of 3D coordinates [x, y, z] for all atoms
            radii (list): List of atomic radii corresponding to each vertex
            vertex_indices (list): Indices of surface atoms to process (subset of vertices)
            probe (float): Radius of the probe sphere used to roll over the surface

        Note:
            - Uses 960 sphere points for surface sampling (generated by asa.generate_sphere_points)
            - Employs spatial hashing (SpaceHash) for efficient neighbor finding
            - Only processes atoms identified as surface atoms (vertex_indices)
            - Test points are placed at distance (atom_radius + probe) from atom centers
        """

        sphere_points = asa.generate_sphere_points(960)
        spacehash = SpaceHash(vertices)
        test_point = [0.0, 0.0, 0.0]

        for i_vertex in tqdm.tqdm(vertex_indices):
            neighbor_indices = spacehash.find_connected_vertex_indices(
                radii, probe, i_vertex
            )
            n_neighbor = len(neighbor_indices)
            i_neighbor_start = 0
            radius_i = probe + radii[i_vertex]
            vertex_i = spacehash.vertices[i_vertex]
            for sphere_point in sphere_points:
                is_point_accessible = True

                test_point[0] = sphere_point[0] * radius_i + vertex_i[0]
                test_point[1] = sphere_point[1] * radius_i + vertex_i[1]
                test_point[2] = sphere_point[2] * radius_i + vertex_i[2]

                for i_neighbor in asa.reorder_range(n_neighbor, i_neighbor_start):
                    j_atom = neighbor_indices[i_neighbor]
                    radius_j = radii[j_atom] + probe
                    vertex_j = spacehash.vertices[j_atom]
                    if pos_distance_sq(test_point, vertex_j) < radius_j * radius_j:
                        i_neighbor_start = i_neighbor
                        is_point_accessible = False
                        break
                if is_point_accessible:
                    self.exclude_sphere(test_point, probe)

    def make_soup(self, res_type, atom_type):
        soup = Soup()
        soup.push_structure_id("HOLLOW")

        element = ""
        for c in atom_type[:2]:
            if not c.isdigit() and c != " ":
                element += c

        i_res = 1
        for i in range(self.n):
            for j in range(self.n):
                for k in range(self.n):
                    if not (self.is_excluded_or_drilled(i, j, k)):
                        pos = self.pos(i, j, k)
                        soup.add_atom(
                            x=pos.x,
                            y=pos.y,
                            z=pos.z,
                            bfactor=0.0,
                            alt="",
                            atom_type=atom_type,
                            elem=element,
                            res_type=res_type,
                            res_num=i_res,
                            ins_code="",
                            chain="A",
                        )
                        i_res += 1
        return soup


def calculate_average_bfactor(grid_chain, protein_atoms, bfactor_probe):
    max_bfactor = 0.0
    for atom in protein_atoms:
        if atom.bfactor > max_bfactor:
            max_bfactor = atom.bfactor
    for grid_atom in grid_chain.atoms():
        bfactors = []
        for protein_atom in protein_atoms:
            if protein_atom.element != "H":
                radius = bfactor_probe
                dist = vector3d.pos_distance(protein_atom.pos, grid_atom.pos)
                if dist < radius:
                    bfactors.append(protein_atom.bfactor)
        n_bfactor = len(bfactors)
        if n_bfactor == 0:
            grid_atom.bfactor = max_bfactor
        else:
            grid_atom.bfactor = sum(bfactors) / float(n_bfactor)


def get_sphere_constraint_fn(center, radius):
    return lambda pos: vector3d.pos_distance(center, pos) <= radius


def get_cylinder_constraint_fn(center1, center2, radius):
    axis12 = center2 - center1

    def cylinder_constraint_fn(pos):
        pos1 = pos - center1
        if vector3d.dot(pos1, axis12) < 0:
            return False
        pos1_perp = pos1.perpendicular_vec(axis12)
        if pos1_perp.mag() > radius:
            return False
        pos2 = pos - center2
        if vector3d.dot(pos2, axis12) > 0:
            return False
        return True

    return cylinder_constraint_fn


def get_constraint(soup, atom_indices, constraint_file, grid_spacing):
    # setup constraints and grid size in width
    constraint_fn = None
    inner_constraint_fn = None
    is_calculate_asa_shell = True

    if not constraint_file:
        center = soup.get_center(atom_indices)
        extent = soup.get_extent_from_center(center, atom_indices)
    else:
        print("Loading constraints from %s" % constraint_file)
        constraints = read_parameters(constraint_file)

        if not constraints.remove_asa_shell:
            is_calculate_asa_shell = False

        if constraints.type == "sphere":
            atom_proxy1 = soup.find_atom_in_soup(
                constraints.chain1, constraints.res_num1, constraints.atom1
            )
            radius = constraints.radius
            constraint_fn = get_sphere_constraint_fn(atom_proxy1.pos, radius)
            center = atom_proxy1.pos
            extent = 2.0 * constraints.radius + 2.0 * grid_spacing
            radius = radius - 3.0 * grid_spacing
            inner_constraint_fn = get_sphere_constraint_fn(atom_proxy1.pos, radius)

        elif constraints.type == "cylinder":
            atom_proxy1 = soup.find_atom_in_soup(
                constraints.chain1, constraints.res_num1, constraints.atom1
            )
            atom_proxy2 = soup.find_atom_in_soup(
                constraints.chain2, constraints.res_num2, constraints.atom2
            )
            axis12 = atom_proxy2.pos - atom_proxy1.pos

            offset1 = -axis12.normalized_vec()
            offset1.scale(constraints.axis_offset1)
            center1 = atom_proxy1.pos + offset1

            offset2 = axis12.normalized_vec()
            offset2.scale(constraints.axis_offset2)
            center2 = atom_proxy2.pos + offset2

            center = center1 + center2
            center.scale(0.5)
            radius = constraints.radius
            constraint_fn = get_cylinder_constraint_fn(center1, center2, radius)

            half_length = vector3d.pos_distance(center, center1)
            extent = 2.0 * grid_spacing + 2.0 * math.sqrt(
                half_length * half_length + constraints.radius * constraints.radius
            )
            border_length = 3.0 * grid_spacing
            center1 = center1 + axis12.normalized_vec().scaled_vec(border_length)
            center2 = center2 - axis12.normalized_vec().scaled_vec(border_length)
            radius = radius - border_length
            inner_constraint_fn = get_cylinder_constraint_fn(center1, center2, radius)

        else:
            raise ValueError("Don't understand constraint type")

    return center, extent, constraint_fn, inner_constraint_fn, is_calculate_asa_shell


def calc_average_bfactor_soup(grid_soup, soup, bfactor_probe):
    """Calculate average B-factors for grid atoms using soup"""
    protein_atom_proxy = soup.get_atom_proxy()
    max_bfactor = 0.0
    for i_protein_atom in range(soup.get_atom_count()):
        protein_atom_proxy.load(i_protein_atom)
        if protein_atom_proxy.bfactor > max_bfactor:
            max_bfactor = protein_atom_proxy.bfactor

    grid_atom_proxy = grid_soup.get_atom_proxy()
    n_grid_atom = grid_soup.get_atom_count()
    for i_grid_atom in tqdm.trange(n_grid_atom):
        grid_atom_proxy.load(i_grid_atom)
        grid_atom_pos = grid_atom_proxy.pos

        bfactors = []
        for i_protein_atom in range(soup.get_atom_count()):
            protein_atom_proxy.load(i_protein_atom)
            if protein_atom_proxy.elem != "H":
                dist = vector3d.pos_distance(protein_atom_proxy.pos, grid_atom_pos)
                if dist < bfactor_probe:
                    bfactors.append(protein_atom_proxy.bfactor)

        n_bfactor = len(bfactors)
        if n_bfactor == 0:
            grid_atom_proxy.bfactor = max_bfactor
        else:
            grid_atom_proxy.bfactor = sum(bfactors) / float(n_bfactor)


def make_hollow_spheres(
    input_file,
    output_file="",
    grid_spacing=defaults.grid_spacing,
    interior_probe=defaults.interior_probe,
    is_skip_waters=defaults.is_skip_waters,
    surface_probe=defaults.surface_probe,
    constraint_file="",
    bfactor_probe=defaults.bfactor_probe,
):
    soup = load_soup(input_file, scrub=True)
    print(
        f"Loaded {soup.get_atom_count()} atoms in {soup.get_residue_count()} residues from `{input_file}`"
    )

    atom_indices = soup.get_atom_indices(skip_waters=is_skip_waters)
    center, extent, constraint_fn, inner_constraint_fn, is_calculate_asa_shell = (
        get_constraint(soup, atom_indices, constraint_file, grid_spacing)
    )
    grid = HollowGrid(grid_spacing, extent, center)
    print(
        f"Grid: {grid.n} x {grid.n} x {grid.n}  Spacing: {grid_spacing:.2f}  Extent: {extent:.2f}"
    )

    print(f"Excluding protein bulk with {interior_probe:.1f} Å probe")
    vertices = []
    radii = []
    atom_proxy = soup.get_atom_proxy()
    for i_atom in atom_indices:
        atom_proxy.load(i_atom)
        vertices.append(atom_proxy.pos.tuple())
        radii.append(atom_proxy.radius)
    grid.exclude_vertices(vertices, radii, interior_probe)

    if constraint_file:
        print("Excluding exterior of constraint")
        grid.exclude_points_in_constraint(constraint_fn)

    if is_calculate_asa_shell:
        # Roll large ball over surface residues, then drill in from the edge
        print("Identifying surface atoms using ASA with 1.4 Å probe")
        atom_asas = asa.calculate_asa_from_soup(soup, 1.4, atom_indices=atom_indices)
        surface_vertex_indices = [i for i, a in enumerate(atom_asas) if a >= 9]

        print(f"Excluding surface shell with {surface_probe:.1f} Å probe")
        grid.exclude_surface(vertices, radii, surface_vertex_indices, surface_probe)

        print("Expansion of exclusion to edge")
        grid.exclude_edge_to_interior()

    print("Excluding encased grid points")
    hole_size = int(1.5 * 1.4 / grid_spacing)
    grid.exclude_surrounded(hole_size)

    # Make hollow spheres from grid-points
    grid_soup = grid.make_soup(defaults.res_type, defaults.atom_type)

    if bfactor_probe:
        print("Averaging nearby protein b-factors for each hollow atom")
        calc_average_bfactor_soup(grid_soup, soup, bfactor_probe)

    if constraint_file:
        print("Setting occupancy to 0 if outside constraint")
        atom_proxy = grid_soup.get_atom_proxy()
        for i_atom in range(grid_soup.get_atom_count()):
            pos = atom_proxy.load(i_atom).pos
            atom_proxy.occupancy = 1.0 if inner_constraint_fn(pos) else 0.0

    if not output_file:
        output_file = add_suffix_to_basename(input_file, "-hollow")
    print("Saving hollow spheres to", output_file)
    write_soup(grid_soup, output_file)


def main():
    @click.command()
    @click.version_option()
    @click.argument("input-file", type=click.Path(exists=True))
    @click.option(
        "-g",
        "--grid-spacing",
        type=float,
        default=defaults.grid_spacing,
        help=f"Grid spacing (default {defaults.grid_spacing:.1f}; 0.2 for final resolution) Å",
    )
    @click.option(
        "-c",
        "--constraint-file",
        type=click.Path(exists=True),
        default=None,
        help="Config file for grid constraints",
    )
    @click.option(
        "-o",
        "--output-file",
        type=click.Path(),
        default="",
        help="Output PDB file for hollow spheres (default: auto-generated)",
    )
    @click.option(
        "-p",
        "--interior-probe",
        type=float,
        default=defaults.interior_probe,
        help=f"Radius of ball to explore cavities (default {defaults.interior_probe:.1f} Å = 95% x radius of output atom type suggested)",
    )
    @click.option(
        "-s",
        "--surface-probe",
        type=float,
        default=defaults.surface_probe,
        help=f"Radius of probe to roll over surface used to define depressions (default {defaults.surface_probe:.2f} angstroms)",
    )
    @click.option(
        "-w",
        "--include-waters",
        is_flag=True,
        default=not defaults.is_skip_waters,
        help="Include water molecules for analysis (default: false)",
    )
    @click.option(
        "-b",
        "--bfactor-probe",
        type=float,
        default=defaults.bfactor_probe,
        help=f"Radius around a grid point, in which the b-factors of heavy atoms are averaged (0.0=off; suggested=4.0; default={defaults.bfactor_probe:.2f})",
    )
    def cli(
        input_file,
        grid_spacing,
        constraint_file,
        output_file,
        interior_probe,
        surface_probe,
        include_waters,
        bfactor_probe,
    ):
        """
        Hollow (c) 2025 Bosco Ho & Franz Gruswitz.

        Generate hollow spheres to fill voids, pockets, clefts and channels in protein structures.

        Creates a PDB file with fake atoms that represent the hollow spaces within the protein
        structure. This is useful for visualizing the surface area of cavities, channels and
        binding sites in Pymol and similar viewers.

        Examples:

            hollow protein.pdb

            hollow protein.pdb --output cavities.pdb

            hollow protein.pdb --grid-spacing 0.3 --interior-probe 1.2
        """
        make_hollow_spheres(
            input_file,
            output_file,
            grid_spacing,
            interior_probe,
            not include_waters,  # Note: inverted logic for is_skip_waters
            surface_probe,
            constraint_file,
            bfactor_probe,
        )

    if len(sys.argv) == 1:
        cli(["--help"])
    else:
        cli()


if __name__ == "__main__":
    main()
