


  # pdbstruct

`pdbstruct` provides a set of tools for fast protein analysis in python. Current tools:

1. `pdbstruct hollow` or `hollow` - generates hollow atoms to facilitate hi-res display of voids, pockets and channels.
2. `pdbstruct volume` - volume calculator and generator of fake atoms that fill the space.
3. `pdstruct asa` - ASA calculator and saves a PDB/CIF file with ASA in the bfactor column.

This was formerly known as [Hollow](https://github.com/boscoh/hollow) but significant improvements resulted in a more general package,
which is capable of forming a foundation for efficient protein analysis:

- modern python3 packaging
- mmCIF parsers and writers
- memory efficient representation of protein
- spatial hashing for fast pair-wise search


## Quick install

1. If you have [uv](https://docs.astral.sh/uv/) installed, then for a global install:

   `>> uv tool install pdbstruct@latest`

2. Another alternative is to use [pipx](https://github.com/pypa/pipx)

   `>> pipx install pdbstruct`

3. If you have a [venv](https://docs.python.org/3/library/venv.html) python environment setup, then:

   `>> pip install pdbstruct`

## Hollow

Hollow was originally developed by Bosco Ho and Franz Gruswitz to solve the problem of displaying protein channels in high resolution. 

There are two effective modes to run Hollow:

1. an automated mode that explicitly deduces the molecular surface
2. a constrained mode that works in a pre-specified volume.

To see the options, type in the command line: 
  
    >> pdbstruct hollow --help

Or for backwards compatibility:

    >> hollow
  
Please refer to the [website](https://boscoh.github.io/pdbstruct/) for illustrated tutorials:

  1) The interior [pathways of myoglobin](https://boscoh.github.io/pdbstruct/myoglobin.html)
  2) The [channel surface]( https://boscoh.github.io/pdbstruct/channel.html) of a phosphate proin:
      

### Automated Surface Detection Mode

In the automated (also the default) mode, a grid is constructed over the entire protein. To analyze large proteins, due to the large size of the resultant grid, use only coarse grid spacings of 1.0 angstroms or 0.5 angstroms. 

The automated mode requires several intermediate calculations. The Accessible Surface Area is first calculated to define surface atoms. Then we sweep a large ball (of say 8.0 angstroms) over all the surface atoms to remove grid points near the surface of the protein. The remaining grid points found under the ball define depressions. This is the costliest part of the calculation.



### Constrained Mode

In the constrained mode, either a sphere or a cylinder is specified in a separate file. Grid points are constructed within these constraints, and thus, there is no need for the costly calculation that sweeps over the entire surface. 

The point of having a constrained mode is that once you identify an interior region of interest from a previous calculation in the Default Mode using a coarse grid, you can identify a smaller region of interest and then set up a very fine grid (eg 0.25 angstroms) around this region.

We also label the occupancy with 1.0 if the sphere is within the accessible surface of the protein, otherwise we label the occupancy with 0.0 if the sphere is outside the accessible surface of the protein. This allows the visualization of funnel-like surfaces for channel proteins.

Here's a sample spherical constraint file. The format is a Python dictionary:

```python
{
  'type': 'sphere',          # 'cylinder' or 'sphere'
  'remove_asa_shell': True,  # True or False
  'radius': 13.0,      
  'chain1': 'A',             
  'res_num1': 107,     
  'atom1': 'CD1',      
}
```

The 'type' refers to whether the constraint is in the shape of a sphere or a cylinder. 

The boolean value 'remove_asa_shell' tells the program whether we want to do the outer shell calculation. If this is True then the fake atoms outside the Accessible Surface will be removed. If False then the fake atoms will follow the shape of the sphere (or cylinder) as it projects out of the exterior of the protein.

The atom around which the spherical constraint is centered is denoted by 'chain1', 'res_num1' and 'atom1'. 

Here's a sample cylinder constraint file:

```python
{
  'type': 'cylinder',          # 'cylinder' or 'sphere'
  'remove_asa_shell': False,   # True or False
  'radius': 10.0,      
  'chain1': 'F',              
  'res_num1': 57,             
  'atom1': 'CE1',             
  'axis_offset1': 0,          
                              
  'chain2': 'G',              
  'res_num2': 185,            
  'atom2': 'CA',              
  'axis_offset2': -3,
}
```

For the cylinder constraint, 'chain1' 'res1' and 'atom1' refers to the center of the one end of the cylinder, whilst 'chain2' 'res2' 'atom2' refers to the center at the other end of the cylinder. The length of the cylinder is inferred from the inter-atomic distance between these two atoms. The 'axis_offset1' and 'axis_offset2' allows the length of the cylinder to be adjusted along the cylinder aixs.



### Default options
  
Default values for various parameters are stored in hollow.txt. If you use the program a lot, you might want to fine tune these options.
            
```python
{
  'grid_spacing': 1.0,
  'interior_probe': 1.444,
  'surface_probe': 8.0,
  'bfactor_probe': 0.0,
  'res_type': 'HOH',
  'atom_type': 'O',
  'atom_field': 'ATOM'
}
```

The 'grid_spacing' determines how detailed the resultant fake atoms will be in the output PDB file. Since the program is written in standard Python, the object that holds the information about the grid runs across memory constraints. That is why it is suggested that a medium resolution of 1.0 angstrom is used in a first spacing. Another problem is that a finer grid will generate a lot of fake atoms. This runs into the problem of displaying the atoms in the protein viewer. 

Fortunately, we have found that the structures at 1.0 and 0.5 angstroms are adequate for identifying most of the regions of interest. And a fine grid of 0.2 angstrom has more than enough detail, when used with constraints.

The definition of the accessible surface involves rolling a surface probe over the atomic surface defined by the van der Waals radius of the atoms. The 'interior probe' defines the size of the surface probe.

In the automatic analysis mode, we want to generate fake atoms in clefts. We start the calculation with a cubic grid of fake atoms and systematically eliminate fake atoms that like within the protein structure. However, we need to also eliminate fake atoms outside the accessible surface of the protein as well, whilst allowing fake atoms to fill in surface clefts. To do that, we need to define a wrapping surface around the protein structure, that wraps around the accessible surface but contains spaces corresponding to surface clefts. We do this by defining a large exterior surface probe through 'surface_probe' (default is 8.0 angstroms).

Three other options are also given, and these relate to the chemistry of the fake atoms in the resultant PDB file. There is also the choice of whether the atoms are written as 'ATOM' or 'HETATM' in the PDB file. The default is 'ATOM' even though water is typically written out as 'HETATM' in order to exploit PyMol. In PyMol, molecular surfaces are only generated for polymers, labeled with 'ATOM'. You can't generate molecular surfaces for 'HETATM' entries.

  ### B-factors

We also calculate appropriate B-factors of every fake atom, by averaging over the heavy protein atoms around each fake atom. This is controlled by the command-line option 'BFACTOR_PROBE'.


  ### Works with PyMol

We developed this program with output designed to be easily viewed and manipulated with PyMol, an open-source protein viewer. By default, the hollow spheres are stored with the "ATOM" field as water oxygen molecules. Pymol can draw the molecular surface of overlapping fake water molecules as it interprets "ATOM" as if the atoms belong to a pseudo polymer.

  ## Change log

- Version 2.0 (Jun 2025). Python 3. pypi. mmCif. Memory effient
    representation of protein. Spatial hashing to speed pair-wise
    search. Removed idle functions.
- Version 1.3 (May 2020). Python 3/2 compatible.</li>
- Version 1.2 (Aug 2011). Changed exceptions to work with Python 2.7
    (thanks Joshua Adelman)
- Version 1.1 (Feb 2009). Faster initialization of grid. Works in the
    IDLE Python interpreter.
