#!/usr/bin/env python
# coding: utf-8

import sys
from typing import Optional

import click

from .asa import calc_asa
from .volume import calc_volume
from .hollow import make_hollow_spheres


@click.group()
@click.version_option()
def cli():
    """PDBStruct - Tools for analyzing protein structures."""
    pass


@cli.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.option(
    "--spacing",
    "-s",
    default=0.5,
    type=float,
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
def volume(
    input_file: str, spacing: float, chain: Optional[str], residue: Optional[int]
):
    """
    Calculate the volume of atoms in a PDB/CIF file using grid-based method.

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
        click.echo("Error: --residue option requires --chain to be specified", err=True)
        sys.exit(1)

    if spacing <= 0:
        click.echo("Error: Grid spacing must be positive", err=True)
        sys.exit(1)

    try:
        calc_volume(input_file, spacing, chain, residue)
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.option(
    "--n-sphere",
    "-n",
    default=960,
    type=int,
    help="Number of points on sphere for calculation (default: 960, more = accurate but slower)",
)
def asa(input_file: str, n_sphere: int):
    """
    Calculate the Accessible Surface Area (ASA) of atoms in a PDB/CIF file.

    Uses the dot density technique with a spherical probe (radius 1.4 Å)
    to calculate the solvent accessible surface area. The algorithm places
    points on a sphere around each atom and tests if they are accessible
    (not buried by neighboring atoms).

    The results are written to a new PDB file with ASA values in the B-factor column.

    Examples:

        pdbstruct asa protein.pdb

        pdbstruct asa protein.pdb --n-sphere 1920
    """
    if n_sphere <= 0:
        click.echo("Error: Number of sphere points must be positive", err=True)
        sys.exit(1)

    try:
        calc_asa(input_file, n_sphere)
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.option(
    "--output",
    "-o",
    type=str,
    help="Output PDB file for hollow spheres (default: auto-generated)",
)
@click.option(
    "--grid-spacing",
    "-g",
    default=0.5,
    type=float,
    help="Grid spacing in Angstroms (default: 0.5)",
)
@click.option(
    "--surface-probe",
    "-s",
    default=1.4,
    type=float,
    help="Radius of surface probe in Angstroms (default: 1.4)",
)
@click.option(
    "--skip-waters/--include-waters",
    default=True,
    help="Skip water molecules in calculation (default: skip)",
)
@click.option(
    "--interior-probe",
    "-p",
    default=1.0,
    type=float,
    help="Radius of interior probe in Angstroms (default: 1.0)",
)
@click.option(
    "--constraint-file",
    "-c",
    type=click.Path(exists=True),
    help="Constraint file for grid constraints",
)
@click.option(
    "--min-volume",
    "-m",
    default=0.0,
    type=float,
    help="Minimum volume threshold (default: 0.0)",
)
def hollow(
    input_file: str,
    output: Optional[str],
    grid_spacing: float,
    surface_probe: float,
    skip_waters: bool,
    interior_probe: float,
    constraint_file: Optional[str],
    min_volume: float,
):
    """
    Generate hollow spheres to fill voids, pockets, clefts and channels in protein structures.

    Creates a PDB file with fake atoms that represent the hollow spaces within
    the protein structure. This is useful for visualizing cavities and binding sites.

    Examples:

        pdbstruct hollow protein.pdb

        pdbstruct hollow protein.pdb --output cavities.pdb

        pdbstruct hollow protein.pdb --grid-spacing 0.3 --interior-probe 1.2
    """
    if grid_spacing <= 0:
        click.echo("Error: Grid spacing must be positive", err=True)
        sys.exit(1)

    if surface_probe <= 0:
        click.echo("Error: Surface probe radius must be positive", err=True)
        sys.exit(1)

    if interior_probe <= 0:
        click.echo("Error: Interior probe radius must be positive", err=True)
        sys.exit(1)

    if min_volume < 0:
        click.echo("Error: Minimum volume must be non-negative", err=True)
        sys.exit(1)

    # Generate default output filename if not provided
    if output is None:
        base_name = input_file.replace(".pdb", "").replace(".cif", "")
        output = f"{base_name}-hollow.pdb"

    try:
        make_hollow_spheres(
            input_file,
            output,
            grid_spacing,
            surface_probe,
            skip_waters,
            interior_probe,
            constraint_file or "",
            min_volume,
        )
        click.echo(f"Hollow spheres written to {output}")
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


if __name__ == "__main__":
    cli()
