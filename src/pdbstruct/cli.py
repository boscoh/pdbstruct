#!/usr/bin/env python
# coding: utf-8

import logging
from typing import Optional

import click

from .asa import calc_asa
from .hollow import make_hollow_spheres
from .util import click_validate_positive, config, init_console_logging
from .volume import calc_volume

logger = logging.getLogger(__name__)


@click.group()
@click.version_option()
def cli():
    """
    PDBStruct - Tools for analyzing protein structures.

    (c) 2025 Bosco Ho & Franz Gruswitz.

    """
    init_console_logging()


@cli.command(no_args_is_help=True)
@click.argument("input_file", type=click.Path(exists=True))
@click.option(
    "--spacing",
    "-s",
    default=0.5,
    type=float,
    callback=click_validate_positive,
    help="Grid spacing in Angstroms (default: 0.5, smaller = more accurate but slower)",
)
@click.option(
    "--chain", "-c", type=str, help="Calculate volume for specific chain only"
)
@click.option(
    "--residue",
    "-r",
    type=int,
    help="Calculate volume for specific residue number in chain (requires --chain)",
)
@click.option(
    "--include-waters",
    default=False,
    is_flag=True,
    help="Include water molecules in calculation (default: False)",
)
def volume(
    input_file: str,
    spacing: float,
    chain: Optional[str],
    residue: Optional[int],
    include_waters: bool,
):
    """
    Calculate the volume of atoms.

    The algorithm uses a 3D grid to discretize space and marks grid points
    that fall within atomic spheres as occupied. The volume is calculated
    as the number of occupied grid points multiplied by the grid spacing cubed.

    Examples:

        pdbstruct volume protein.pdb

        pdbstruct volume protein.pdb --spacing 0.3

        pdbstruct volume protein.pdb --chain A

        pdbstruct volume protein.pdb --chain A --residue 100
    """
    if residue is not None and chain is None:
        raise click.UsageError(
            "Error: --residue option requires --chain to be specified"
        )
    calc_volume(input_file, spacing, chain, residue, not include_waters)


@cli.command(no_args_is_help=True)
@click.argument("input_file", type=click.Path(exists=True))
@click.option(
    "--n-sphere",
    "-n",
    default=960,
    type=int,
    callback=click_validate_positive,
    help="Number of points on sphere for calculation (default: 960, more = accurate but slower)",
)
@click.option(
    "--include-waters",
    default=False,
    is_flag=True,
    help="Include water molecules in calculation (default: False)",
)
def asa(input_file: str, n_sphere: int, include_waters: bool):
    """
    Calculate the accessible-surface-area (ASA) of atoms.

    Uses the dot density technique with a spherical probe (radius 1.4 Å)
    to calculate the solvent accessible surface area. The algorithm places
    points on a sphere around each atom and tests if they are accessible
    (not buried by neighboring atoms).

    The results are written to a new PDB file with ASA values in the B-factor column.

    Examples:

        pdbstruct asa protein.pdb

        pdbstruct asa protein.pdb --n-sphere 1920
    """
    calc_asa(input_file, n_sphere, not include_waters)


@cli.command(no_args_is_help=True)
@click.version_option()
@click.argument("input-file", type=click.Path(exists=True))
@click.option(
    "-g",
    "--grid-spacing",
    type=float,
    default=config.grid_spacing,
    callback=click_validate_positive,
    help=f"Grid spacing (default {config.grid_spacing:.1f}; 0.2 for final resolution) Å",
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
    default=config.interior_probe,
    callback=click_validate_positive,
    help=f"Radius of ball to explore cavities (default {config.interior_probe:.1f} Å = 95% x radius of output atom type suggested)",
)
@click.option(
    "-s",
    "--surface-probe",
    type=float,
    default=config.surface_probe,
    callback=click_validate_positive,
    help=f"Radius of probe to roll over surface used to define depressions (default {config.surface_probe:.2f} angstroms)",
)
@click.option(
    "-w",
    "--include-waters",
    is_flag=True,
    default=not config.is_skip_waters,
    help="Include water molecules for analysis (default: false)",
)
@click.option(
    "-b",
    "--bfactor-probe",
    type=float,
    default=config.bfactor_probe,
    callback=click_validate_positive,
    help=f"Radius around a grid point, in which the b-factors of heavy atoms are averaged (0.0=off; suggested=4.0; default={config.bfactor_probe:.2f})",
)
@click.option(
    "--pymol",
    type=bool,
    default=False,
    is_flag=True,
    help="Show hollow spheres in pymol (if found)",
)
def hollow(
    input_file: str,
    output_file: Optional[str],
    grid_spacing: float,
    surface_probe: float,
    include_waters: bool,
    interior_probe: float,
    constraint_file: Optional[str],
    bfactor_probe: float,
    pymol: bool,
):
    """
    Generate spheres to fill voids, pockets and channels.

    Creates a PDB file with fake atoms that represent the hollow spaces within the protein
    structure. This is useful for visualizing the surface area of cavities, channels and
    binding sites in Pymol and similar viewers.

    Examples:

        pdbstruct hollow protein.pdb

        pdbstruct hollow protein.pdb --output_file cavities.pdb

        pdbstruct hollow protein.pdb --grid-spacing 0.3 --interior-probe 1.2
    """
    make_hollow_spheres(
        input_file=input_file,
        output_file=output_file,
        grid_spacing=grid_spacing,
        interior_probe=interior_probe,
        is_skip_waters=not include_waters,
        surface_probe=surface_probe,
        constraint_file=constraint_file or "",
        bfactor_probe=bfactor_probe,
        show_pymol=pymol,
    )


if __name__ == "__main__":
    cli()
