import array
import math
from typing import Generator, List, Tuple

from .vector3d import pos_distance_sq as pos_distance_sq


class SpaceHash(object):
    """
    A spatial hash data structure for efficient 3D proximity queries on vertex collections.

    SpaceHash partitions 3D space into a regular grid of voxels and maps vertices to hash buckets
    based on their spatial coordinates. This enables efficient nearest neighbor searches and
    proximity queries by limiting searches to spatially adjacent voxels, reducing complexity
    from O(n²) to approximately O(n) for uniformly distributed data.

    The hash structure is particularly useful for molecular modeling applications where you need
    to frequently find atoms within a certain distance of each other, such as for contact
    detection, surface calculations, or bond identification.

    Attributes:
        vertices: List of 3D vertex coordinates [[x, y, z], ...]
        div (float): Voxel size for spatial partitioning (default 5.3)
        inv_div (float): Inverse of div for efficient division operations
        padding (float): Extra padding added to bounding box (default 0.05)
        minima (List[float]): Minimum coordinates in each dimension [x_min, y_min, z_min]
        maxima (List[float]): Maximum coordinates in each dimension [x_max, y_max, z_max]
        spans (List[float]): Span in each dimension (maxima - minima)
        sizes (List[int]): Number of voxels in each dimension
        vertex_indices_by_hash (dict): Maps hash values to arrays of vertex indices
        voxels (List[List[int]]): Voxel coordinates for each vertex

    Args:
        vertex_indices: List of 3D coordinates representing vertices/atoms
        div (float, optional): Size of each voxel in the spatial hash grid. Smaller values
                              create finer resolution but use more memory. Defaults to 5.3.
        padding (float, optional): Extra padding added around the bounding box to handle
                                  edge cases. Defaults to 0.05.

    Example:
        >>> vertices = [[0, 0, 0], [1, 1, 1], [5, 5, 5], [1.2, 0.8, 1.1]]
        >>> space_hash = SpaceHash(vertices, div=2.0)
        >>>
        >>> # Find all close pairs
        >>> for i, j in space_hash.close_pairs():
        ...     print(f"Vertices {i} and {j} are spatially close")
        >>>
        >>> # Find atoms connected to atom 0 within probe distance
        >>> radii = [1.0, 1.2, 1.5, 1.1]  # Atomic radii
        >>> connected = space_hash.find_connected_vertex_indices(radii, probe=1.4, i_vertex=0)
        >>> print(f"Atoms connected to atom 0: {connected}")

    Note:
        The choice of div parameter affects both performance and accuracy. Too large values
        may miss nearby vertices in adjacent voxels, while too small values increase memory
        usage and reduce the efficiency gains of spatial hashing.
    """

    def __init__(self, vertex_indices, div=5.3, padding=0.05):
        self.vertices = vertex_indices
        self.div = div
        self.inv_div = 1.0 / self.div
        self.padding = padding

        zero3 = lambda: [0.0] * 3
        self.minima = zero3()
        self.maxima = zero3()
        self.spans = zero3()
        self.sizes = zero3()

        for i in range(3):
            self.minima[i] = min([v[i] for v in self.vertices])
            self.maxima[i] = max([v[i] for v in self.vertices])
            self.minima[i] -= self.padding
            self.maxima[i] += self.padding
            self.spans[i] = self.maxima[i] - self.minima[i]
            self.sizes[i] = int(math.ceil(self.spans[i] * self.inv_div))

        self.size1_size2 = self.sizes[1] * self.sizes[2]
        self.size2 = self.sizes[2]

        self.vertex_indices_by_hash = {}
        self.voxels = []
        for i_vertex, vertex in enumerate(self.vertices):
            voxel = self.vertex_to_voxel(vertex)
            self.voxels.append(voxel)
            hash = self.voxel_to_hash(voxel)
            vertex_indices = self.vertex_indices_by_hash.setdefault(
                hash, array.array("L")
            )
            vertex_indices.append(i_vertex)

    def vertex_to_voxel(self, v):
        return [int((v[i] - self.minima[i]) * self.inv_div) for i in range(3)]

    def voxel_to_hash(self, s):
        return s[0] * self.size1_size2 + s[1] * self.size2 + s[2]

    def get_neighbour_voxels(self, voxel):
        def neighbourhood_in_dim(space, i_dim):
            i = max(0, space[i_dim] - 1)
            j = min(self.sizes[i_dim], space[i_dim] + 2)
            return range(i, j)

        for s0 in neighbourhood_in_dim(voxel, 0):
            for s1 in neighbourhood_in_dim(voxel, 1):
                for s2 in neighbourhood_in_dim(voxel, 2):
                    yield [s0, s1, s2]

    def get_vertex_indices_in_voxel(self, voxel):
        hash = self.voxel_to_hash(voxel)
        return self.vertex_indices_by_hash.get(hash, [])

    def close_pairs(self) -> Generator[Tuple[int, int], None, None]:
        """
        Generator that yields pairs of vertex indices that are spatially close to each other.

        This method uses the spatial hash structure to efficiently find pairs of vertices
        that are in the same or neighboring hash cells, avoiding the need to check all
        possible pairs (O(n²) complexity). Only pairs where the first index is less than
        the second index are yielded to avoid duplicates.

        Yields:
            tuple[int, int]: A tuple containing two vertex indices (i_vertex0, i_vertex1)
                            where i_vertex0 < i_vertex1, representing vertices that are
                            spatially close to each other based on the hash cell structure.

        Note:
            This method identifies potentially close pairs based on spatial proximity
            within the hash grid. The actual distance between vertices should be
            calculated separately if exact distance thresholds are needed.

        Example:
            >>> space_hash = SpaceHash(vertices)
            >>> for i, j in space_hash.close_pairs():
            ...     distance = calculate_distance(vertices[i], vertices[j])
            ...     if distance < threshold:
            ...         # Process close pair
        """
        for i_vertex in range(len(self.vertices)):
            voxel_i = self.voxels[i_vertex]
            for voxel_j in self.get_neighbour_voxels(voxel_i):
                for j_vertex in self.get_vertex_indices_in_voxel(voxel_j):
                    if i_vertex < j_vertex:
                        yield i_vertex, j_vertex

    def find_connected_vertex_indices(
        self, radii: List[float], probe: float, i_vertex: int
    ) -> List[int]:
        """
        Returns list of indices of atoms within probe distance to atom k, given the radii of the atoms.

        Args:
            radii (List[float]): List of atomic radii for all vertices/atoms in the structure
            probe (float): Probe radius to add to the search distance (typically solvent probe radius)
            i_vertex (int): Index of the target vertex/atom to find connections for

        Returns:
            List[int]: List of vertex indices that are within the specified distance of the target atom

        Note:
            The method searches for atoms within a distance of (radii[i_vertex] + probe + probe + radii[j_vertex])
            from the target atom. This accounts for the radii of both atoms plus twice the probe radius.

        Example:
            >>> vertices = [[0, 0, 0], [1, 1, 1], [5, 5, 5]]
            >>> space_hash = SpaceHash(vertices)
            >>> atomic_radii = [1.0, 1.2, 1.5]  # Radii for each atom
            >>> probe_radius = 1.4  # Water probe radius
            >>> connected = space_hash.find_connected_vertex_indices(atomic_radii, probe_radius, 0)
            >>> print(f"Atoms connected to atom 0: {connected}")
        """
        voxel_i = self.voxels[i_vertex]
        vertex_indices = list(self.get_vertex_indices_in_voxel(voxel_i))
        for voxel_j in self.get_neighbour_voxels(voxel_i):
            vertex_indices.extend(self.get_vertex_indices_in_voxel(voxel_j))

        result = []

        r0 = radii[i_vertex] + probe + probe
        for j_vertex in vertex_indices:
            if j_vertex == i_vertex:
                continue
            r = r0 + radii[j_vertex]
            if pos_distance_sq(self.vertices[i_vertex], self.vertices[j_vertex]) < r * r:
                result.append(j_vertex)

        return result
