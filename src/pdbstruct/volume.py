#!/usr/bin/env python
# coding: utf-8

import logging
import math
import sys
from typing import List, Tuple

from .bgrid import BoolGrid
from .parse import add_suffix_to_basename, load_soup, write_soup
from .soup import Soup
from .util import tqdm_iter, tqdm_range
from .vector3d import Vector3d

logger = logging.getLogger(__name__)


class VolumeGrid(BoolGrid):
    """3D grid for volume calculations."""

    def __init__(self, grid_spacing: float, width: float, center: Vector3d):
        self.width = float(width)
        half_width = self.width / 2.0
        self.center = center.copy()
        self.spacing = float(grid_spacing)
        self.inv_spacing = 1.0 / self.spacing

        n = 1
        cover = 0
        self.low = Vector3d()

        # Calculate grid size to cover the required width
        while cover < half_width:
            n += 1
            half_n_point = int(n / 2)
            self.low.x = self.center.x - half_n_point * self.spacing
            self.low.y = self.center.y - half_n_point * self.spacing
            self.low.z = self.center.z - half_n_point * self.spacing
            width_1 = abs(self.center.x - self.low.x)
            high_x = self.low.x + n * self.spacing
            width_2 = abs(high_x - self.center.x)
            cover = min(width_1, width_2)

        # Initialize BooleanGrid with calculated n
        super().__init__(n)

        self.actual_width = self.n * self.spacing
        self.total_volume = self.actual_width**3
        self.total_point = self.n**3

        # Pre-calculate coordinate arrays for faster access
        self.x = [self.low.x + i * self.spacing for i in range(self.n)]
        self.y = [self.low.y + j * self.spacing for j in range(self.n)]
        self.z = [self.low.z + k * self.spacing for k in range(self.n)]

    def indices(self, pos: Vector3d) -> Tuple[float, float, float]:
        """Convert 3D position to grid indices (as floats)."""
        return (
            (pos.x - self.low.x) * self.inv_spacing,
            (pos.y - self.low.y) * self.inv_spacing,
            (pos.z - self.low.z) * self.inv_spacing,
        )

    def pos(self, i: int, j: int, k: int) -> Vector3d:
        """Get 3D position from grid indices."""
        return Vector3d(self.x[i], self.y[j], self.z[k])

    def is_grid_point_near_sphere(
        self, i: int, j: int, k: int, pos: Vector3d, r_sq: float
    ) -> bool:
        """Check if grid point (i,j,k) is within sphere at pos with radius squared r_sq."""
        d_x = self.x[i] - pos.x
        d_y = self.y[j] - pos.y
        d_z = self.z[k] - pos.z
        return (d_x * d_x + d_y * d_y + d_z * d_z) < r_sq

    def int_range(self, low_f: float, high_f: float) -> List[int]:
        """Convert float range to integer range with bounds checking."""
        low = max(0, int(math.floor(low_f - 1)))
        high = min(self.n, int(math.ceil(high_f) + 2))
        return list(range(low, high))

    def exclude_sphere(self, pos: Vector3d, r: float):
        """Mark grid points within sphere as excluded (occupied)."""
        low = Vector3d(pos.x - r, pos.y - r, pos.z - r)
        low_i, low_j, low_k = self.indices(low)
        high = Vector3d(pos.x + r, pos.y + r, pos.z + r)
        high_i, high_j, high_k = self.indices(high)
        r_sq = r * r

        for i in self.int_range(low_i, high_i):
            for j in self.int_range(low_j, high_j):
                for k in self.int_range(low_k, high_k):
                    if not self.is_set(i, j, k):  # If not already excluded
                        if self.is_grid_point_near_sphere(i, j, k, pos, r_sq):
                            self.set(i, j, k, True)

    def n_excluded(self) -> int:
        """Count number of excluded (occupied) grid points."""
        return self.count_set()

    def make_soup(self, res_type: str = "HOH", atom_type: str = "O") -> Soup:
        """Create a Soup object with atoms at excluded grid points."""
        soup = Soup()
        soup.push_structure_id("GRID", "Volume Grid Visualization")

        # Determine element from atom type
        element = ""
        for c in atom_type[:2]:
            if not c.isdigit() and c != " ":
                element += c
        if not element:
            element = "O"  # Default to oxygen

        i_res = 1
        for i in tqdm_range(self.n):
            for j in range(self.n):
                for k in range(self.n):
                    if self.is_set(i, j, k):  # If grid point is excluded
                        grid_pos = self.pos(i, j, k)
                        soup.add_atom(
                            x=grid_pos.x,
                            y=grid_pos.y,
                            z=grid_pos.z,
                            bfactor=50.0,  # Default B-factor
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


def calculate_volume_of_atoms(
    soup: Soup,
    atom_indices: List[int],
    grid_spacing: float = 0.5,
    out_fname: str = "",
) -> float:
    """
    Calculate the volume of atoms in a Soup using grid-based method.

    Args:
        soup: Soup object containing atoms
        grid_spacing: Grid spacing in Angstroms (smaller = more accurate but slower)
        atom_indices: List of atom indices to include
        out_fname: Output filename for grid visualization (empty = no output)

    Returns:
        Tuple of (volume in Å³, number of excluded grid points)
    """
    # Calculate center and width
    center = soup.get_center()
    extent = soup.get_extent_from_center(center, atom_indices)

    # Create grid
    grid = VolumeGrid(grid_spacing, extent, center)
    logger.info(
        f"Grid: {grid.n}³;  {grid.n} x {grid_spacing}Å = {grid.actual_width:.1f}Å"
    )

    # Exclude spheres for each atom
    atom_proxy = soup.get_atom_proxy()
    logger.info("Excluding empty space from grid")
    for i_atom in tqdm_iter(atom_indices):
        atom_proxy.load(i_atom)
        grid.exclude_sphere(atom_proxy.pos, atom_proxy.radius)

    logger.info("Counting volume voxels")
    n_voxel = 0
    for i in tqdm_range(grid.n**3):
        if grid.bits[i]:
            n_voxel += 1
    volume_voxel = float(grid_spacing) ** 3
    volume = n_voxel * volume_voxel
    logger.info(f"Volume: {volume:.2f} Å³")

    if out_fname:
        logger.info("Generating volume atoms")
        grid_soup = grid.make_soup("HOH", "O")
        write_soup(grid_soup, out_fname)

    return volume


def calculate_volume_by_residue(
    soup: Soup, i_res: int, grid_spacing: float = 0.5, out_fname: str = ""
) -> float:
    """
    Calculate volume of a specific residue.

    Args:
        soup: Soup object
        i_res: Residue index
        grid_spacing: Grid spacing in Angstroms
        out_fname: Output filename for visualization

    Returns:
        Tuple of (volume in Å³, number of excluded grid points)
    """
    if i_res >= soup.get_residue_count():
        raise ValueError(f"Residue index {i_res} out of range")

    residue_proxy = soup.get_residue_proxy(i_res)
    atom_indices = residue_proxy.get_atom_indices()

    logger.info(
        f"Calculating volume for residue {residue_proxy.res_type} {residue_proxy.res_num} "
        f"chain {residue_proxy.chain} ({len(atom_indices)} atoms)"
    )

    return calculate_volume_of_atoms(soup, atom_indices, grid_spacing, out_fname)


def calculate_volume_by_chain(
    soup: Soup, chain: str, grid_spacing: float = 0.5, out_fname: str = ""
) -> float:
    """
    Calculate volume of all atoms in a specific chain.

    Args:
        soup: Soup object
        chain: Chain identifier
        grid_spacing: Grid spacing in Angstroms
        out_fname: Output filename for visualization

    Returns:
        Tuple of (volume in Å³, number of excluded grid points)
    """
    atom_indices = []
    residue_proxy = soup.get_residue_proxy()

    for i_res in range(soup.get_residue_count()):
        residue_proxy.load(i_res)
        if residue_proxy.chain == chain:
            atom_indices.extend(residue_proxy.get_atom_indices())

    if not atom_indices:
        logger.warning(f"No atoms found in chain {chain}")
        return 0.0

    logger.info(f"Calculating volume for chain {chain} ({len(atom_indices)} atoms)")

    return calculate_volume_of_atoms(soup, atom_indices, grid_spacing, out_fname)


def calc_volume(
    input_file, spacing=0.5, target_chain=None, target_res_num=None, skip_waters=False
) -> float:
    soup = load_soup(input_file, scrub=True)

    if target_chain and target_res_num is not None:
        # Find specific residue in chain
        residue_proxy = soup.get_residue_proxy()
        target_res_index = None

        for i_res in range(soup.get_residue_count()):
            residue_proxy.load(i_res)
            if (
                residue_proxy.chain == target_chain
                and residue_proxy.res_num == target_res_num
            ):
                target_res_index = i_res
                break

        if target_res_index is None:
            logger.error(f"Residue {target_res_num} not found in chain {target_chain}")
            sys.exit(1)

        out_fname = add_suffix_to_basename(
            input_file, f"-chain{target_chain}-res{target_res_num}-volume"
        )
        volume = calculate_volume_by_residue(soup, target_res_index, spacing, out_fname)

    elif target_chain:
        # Calculate volume for entire chain
        out_fname = add_suffix_to_basename(input_file, f"-chain{target_chain}-volume")
        volume = calculate_volume_by_chain(soup, target_chain, spacing, out_fname)

    else:
        # Calculate total volume
        atom_indices = soup.get_atom_indices(skip_waters=skip_waters)
        out_fname = add_suffix_to_basename(input_file, "-volume")
        volume = calculate_volume_of_atoms(soup, atom_indices, spacing, out_fname)

    return volume
