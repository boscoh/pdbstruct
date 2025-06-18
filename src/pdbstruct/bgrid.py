try:
    from bitarray import bitarray
except ImportError:
    raise ImportError(
        "bitarray package is required. Install with: pip install bitarray"
    )


class BoolGrid:
    """
    BoolGrid - A memory-efficient 3D boolean grid using bitarray.

    BoolGrid provides a compact representation of a 3D boolean grid by using
    bitarray storage. Each grid point stores a boolean state (occupied/unoccupied)
    using only 1 bit per point, making it 8x more memory-efficient than the previous
    array-based implementation and highly suitable for large 3D grids used in
    molecular modeling, voxel operations, and spatial data structures.

    The grid uses a flattened 1D bitarray internally with 3D indexing via the formula:
    index = i * n² + j * n + k, where (i,j,k) are the 3D coordinates.

    Common Usage:

    1. Creating a 3D grid for molecular volume calculations:

        # Create a 100x100x100 grid
        grid = BoolGrid(100)

        # Mark some points as occupied
        grid.set(50, 50, 50, True)   # Center point
        grid.set(25, 25, 25, True)   # Another point

        # Check if points are occupied
        if grid.is_set(50, 50, 50):
            print("Center point is occupied")


    2. Using for collision detection or space partitioning:

        # Create grid for spatial partitioning
        grid = BoolGrid(64)

        # Mark occupied regions
        for x in range(10, 20):
            for y in range(10, 20):
                for z in range(10, 20):
                    grid.set(x, y, z, True)

        # Check if region is free
        def is_region_free(start_x, start_y, start_z, size):
            for i in range(start_x, start_x + size):
                for j in range(start_y, start_y + size):
                    for k in range(start_z, start_z + size):
                        if grid.is_set(i, j, k):
                            return False
            return True


    3. Resetting and reusing the grid:

        grid = BoolGrid(50)

        # Use grid for first calculation
        grid.set(10, 10, 10, True)
        grid.set(20, 20, 20, True)

        # Clear and reuse for next calculation
        grid.reset()

        # Grid is now empty and ready for reuse
        assert not grid.is_set(10, 10, 10)


    4. Integration with molecular structures:

        # Grid for protein volume calculation
        spacing = 0.5  # Angstroms
        width = 50.0   # Angstroms
        n = int(width / spacing)  # 100 grid points

        grid = BoolGrid(n)

        # Mark grid points occupied by atoms
        for atom in protein_atoms:
            # Convert atom position to grid coordinates
            i = int((atom.x - min_x) / spacing)
            j = int((atom.y - min_y) / spacing)
            k = int((atom.z - min_z) / spacing)

            # Mark sphere around atom as occupied
            radius_grid = int(atom.radius / spacing)
            for di in range(-radius_grid, radius_grid + 1):
                for dj in range(-radius_grid, radius_grid + 1):
                    for dk in range(-radius_grid, radius_grid + 1):
                        if di*di + dj*dj + dk*dk <= radius_grid*radius_grid:
                            grid.set(i + di, j + dj, k + dk, True)


    Memory Usage (compared to previous array implementation):
    - A 100³ grid uses ~125KB of memory (1,000,000 bits ≈ 125KB vs 1MB previously)
    - A 200³ grid uses ~1MB of memory (8,000,000 bits ≈ 1MB vs 8MB previously)
    - 8x more memory-efficient than the previous byte-based array implementation
    - Much more memory-efficient than using nested lists or numpy boolean arrays

    Performance Notes:
    - Bounds checking is performed on all operations for safety
    - Out-of-bounds access returns False for is_set() and is ignored for set()
    - Direct bitarray access is used internally for maximum speed
    - Reset operation is optimized to clear the entire bitarray efficiently
    - Bitwise operations are available for advanced use cases

    Attributes:
        n (int): Size of each dimension (grid is n×n×n)
        n_sq (int): n², cached for performance
        n_cube (int): n³, total number of grid points
        bits (bitarray): Internal 1D bitarray storing boolean values as bits
    """

    def __init__(self, n: int):
        self.n = n
        self.n_sq = self.n * self.n
        self.n_cube = self.n_sq * self.n
        # Initialize bitarray with all bits set to False (0)
        self.bits = bitarray(self.n_cube)
        self.bits.setall(0)

    def is_set(self, i: int, j: int, k: int) -> bool:
        """Check if grid point is set"""
        if not (0 <= i < self.n and 0 <= j < self.n and 0 <= k < self.n):
            return False
        index = i * self.n_sq + j * self.n + k
        return bool(self.bits[index])

    def set(self, i: int, j: int, k: int, is_state: bool):
        """Set grid point state"""
        if not (0 <= i < self.n and 0 <= j < self.n and 0 <= k < self.n):
            return
        index = i * self.n_sq + j * self.n + k
        self.bits[index] = is_state

    def reset(self):
        """Reset all grid points to unoccupied."""
        self.bits.setall(0)

    def count_set(self) -> int:
        """Count the number of set (True) bits in the grid."""
        return self.bits.count()

    def count_unset(self) -> int:
        """Count the number of unset (False) bits in the grid."""
        return self.n_cube - self.bits.count()

    def flip(self, i: int, j: int, k: int):
        """Flip the state of a grid point."""
        if not (0 <= i < self.n and 0 <= j < self.n and 0 <= k < self.n):
            return
        index = i * self.n_sq + j * self.n + k
        self.bits[index] = not self.bits[index]

    def flip_all(self):
        """Flip all bits in the grid."""
        self.bits.invert()
