
  # pdbstruct

`pdbstruct` - tools for efficient protein analysis in python:

1. `pdbstruct hollow` - generates hollow atoms for hi-res display of voids, pockets and channels.
2. `pdbstruct volume` - volume calculator and generator of fake atoms that fill the space.
3. `pdstruct asa` - calculates and saves atomic ASA to the bfactor column.

This was formerly known as [Hollow](https://github.com/boscoh/hollow) but was renamed because significant improvements in v2.0 means the package can serve as a general module for efficient protein analysis:

- modern packaging and cli
- mmCIF parsers and writers
- memory efficient representation of protein
- spatial hashing for fast close-pair search
- 3d boolean grid


## Quick install

1. If you have [uv](https://docs.astral.sh/uv/) installed, then for a global install:

       >> uv tool install pdbstruct@latest

2. Or you can use ux to run the command in an isolated enviornment:

       >> uvx pdbstruct

4. Another alternative is to use [pipx](https://github.com/pypa/pipx) to install a gloabl cli:

       >> pipx install pdbstruct

5. Or drop into your [venv](https://docs.python.org/3/library/venv.html) python environment:

       >> pip install pdbstruct

  ## Hollow

Hollow was originally developed by Bosco Ho and Franz Gruswitz to solve the problem of displaying protein channels in high resolution. Read more about [Hollow](https://boscoh.github.io/hollow/).

  ## Change log

- Version 2.0 (Jun 2025). Renamed to pdbstruct. Python 3. Pypi. MmCif. Memory effient
    representation of protein. Spatial hashing to speed pair-wise
    search. Removed idle functions.
- Version 1.3 (May 2020). Python 3/2 compatible.</li>
- Version 1.2 (Aug 2011). Changed exceptions to work with Python 2.7
    (thanks Joshua Adelman)
- Version 1.1 (Feb 2009). Faster initialization of grid. Works in the
    IDLE Python interpreter.
