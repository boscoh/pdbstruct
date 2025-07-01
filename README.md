
  # pdbstruct

`pdbstruct` provides a set of tools for fast protein analysis in python:

1. `pdbstruct hollow` - generates hollow atoms to facilitate hi-res display of voids, pockets and channels.
2. `pdbstruct volume` - volume calculator and generator of fake atoms that fill the space.
3. `pdstruct asa` - ASA calculator and saves atomic ASA to the bfactor column.

This was formerly known as [Hollow](https://github.com/boscoh/hollow) but significant improvements resulted in a more general package for efficient protein analysis:

- modern python3 packaging
- mmCIF parsers and writers
- memory efficient representation of protein
- spatial hashing for fast pair-wise search


## Quick install

1. If you have [uv](https://docs.astral.sh/uv/) installed, then for a global install:

       >> uv tool install pdbstruct@latest

2. Another alternative is to use [pipx](https://github.com/pypa/pipx)

       >> pipx install pdbstruct

3. If you have a [venv](https://docs.python.org/3/library/venv.html) python environment setup, then:

       >> pip install pdbstruct

  ## Hollow

Hollow was originally developed by Bosco Ho and Franz Gruswitz to solve the problem of displaying protein channels in high resolution. 

To read more about Hollow, please refer to [here](https://boscoh.github.io/hollow/).

  ## Change log

pdbstruct is a rename of hollow, as the latest version of hollow can serve
as a base for efficient tools for protein analysis.

- Version 2.0 (Jun 2025). Python 3. pypi. mmCif. Memory effient
    representation of protein. Spatial hashing to speed pair-wise
    search. Removed idle functions.
- Version 1.3 (May 2020). Python 3/2 compatible.</li>
- Version 1.2 (Aug 2011). Changed exceptions to work with Python 2.7
    (thanks Joshua Adelman)
- Version 1.1 (Feb 2009). Faster initialization of grid. Works in the
    IDLE Python interpreter.
