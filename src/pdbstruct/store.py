import array
from typing import Callable, List, Optional, Tuple


def get_typed_array(dtype: str, size: int) -> array.array:
    """Get a Python array of the specified type and size"""
    type_mapping = {
        "int8": "b",  # signed char
        "int16": "h",  # signed short
        "int32": "i",  # signed int
        "uint8": "B",  # unsigned char
        "uint16": "H",  # unsigned short
        "uint32": "I",  # unsigned int
        "float32": "f",  # float
    }
    array_type = type_mapping.get(dtype, "f")
    return array.array(array_type, [0] * size)


class Store:
    """
    Store - an array-based space-efficient typed datastore using the Flyweight pattern.

    The Store class provides a columnar data structure where each field is stored as a
    separate typed Python array. This design is memory-efficient for storing large numbers
    of records with the same structure, such as particles, vertices, or database rows.

    Common Usage:

    1. Creating a store for 3D points with colors:
        ```python
        fields = [
            ("x", "float32"), ("y", "float32"), ("z", "float32"),  # coordinates
            ("r", "uint8"), ("g", "uint8"), ("b", "uint8")         # RGB colors
        ]
        store = Store(fields, initial_size=1000)
        ```

    2. Adding data to the store:
        ```python
        # Add a red point at (1.0, 2.0, 3.0)
        index = store.count
        store.increment()  # Grows store if needed
        store.x[index] = 1.0
        store.y[index] = 2.0
        store.z[index] = 3.0
        store.r[index] = 255
        store.g[index] = 0
        store.b[index] = 0
        ```

    3. Sorting by distance from origin:
        ```python
        def compare_by_distance(i, j):
            dist_i = store.x[i]**2 + store.y[i]**2 + store.z[i]**2
            dist_j = store.x[j]**2 + store.y[j]**2 + store.z[j]**2
            return int(dist_i - dist_j)

        store.sort(compare_by_distance)
        ```

    4. Copying data between stores:
        ```python
        # Copy first 10 records from source_store to target_store
        target_store.copy_from(source_store, 0, 0, 10)
        ```

    5. Dynamic field addition:
        ```python
        # Add velocity field after creation
        store.add_field("velocity", "float32")
        ```

    The Store automatically manages memory allocation, growing by 1.5x when full
    (minimum 256 entries). Fields are accessed as attributes (e.g., store.x[0])
    and support all Python array operations on the underlying data.

    Typical use cases:
    - Particle systems (position, velocity, color, age)
    - 3D mesh vertices (position, normal, UV coordinates)
    - Database-like records with mixed data types
    - Time series data with multiple channels
    - Game entities (position, health, score, etc.)
    """

    def __init__(self, fields: List[Tuple[str, str]], size: Optional[int] = None):
        """
        Initialize the Store

        Args:
            fields: list of typed fields in the store (name, type)
            size: initial size of the datastore
        """
        # actual size allocated
        self.capacity = 0

        # size to use, updated by grow_if_full
        self.count = 0
        self._fields = fields

        if isinstance(size, int):
            self._init(size)
        else:
            self._init(0)

    def _init(self, size: int):
        """
        Initialize the store

        Args:
            size: size to initialize
        """
        self.capacity = size
        self.count = 0
        for field in self._fields:
            self._init_field(*field)

    def _init_field(self, name: str, dtype: str):
        """
        Initialize a field

        Args:
            name: field name
            dtype: data type, one of int8, int16, int32, uint8, uint16, uint32, float32
        """
        setattr(self, name, get_typed_array(dtype, self.capacity))

    def add_field(self, name: str, dtype: str):
        """
        Add a field

        Args:
            name: field name
            dtype: data type, one of int8, int16, int32, uint8, uint16, uint32, float32
        """
        self._fields.append((name, dtype))
        self._init_field(name, dtype)

    def resize(self, size: int):
        """
        Resize the store to the new size

        Args:
            size: new size
        """
        self.capacity = round(size or 0)
        self.count = min(self.count, self.capacity)

        for name, dtype in self._fields:
            old_array = getattr(self, name)
            new_array = get_typed_array(dtype, self.capacity)

            # Copy existing data to new array
            copy_length = min(len(old_array), self.capacity)
            for i in range(copy_length):
                new_array[i] = old_array[i]

            setattr(self, name, new_array)

    def grow_if_full(self):
        """
        Resize the store to 1.5 times its current size if full, or to 256 if empty
        """
        if self.count >= self.capacity:
            size = round(self.capacity * 1.5)
            self.resize(max(256, size))

    def increment(self):
        """
        Increase the store by 1, and resize if necessary
        """
        self.count += 1
        self.grow_if_full()

    def copy_from(
        self, other: "Store", this_offset: int, other_offset: int, length: int
    ):
        """
        Copy data from one store to another

        Args:
            other: store to copy from
            this_offset: offset to start copying to
            other_offset: offset to start copying from
            length: number of entries to copy
        """
        for name, _ in self._fields:
            this_field = getattr(self, name)
            other_field = getattr(other, name)

            for j in range(length):
                this_field[this_offset + j] = other_field[other_offset + j]

    def copy_within(self, offset_target: int, offset_source: int, length: int):
        """
        Copy data within this store

        Args:
            offset_target: offset to start copying to
            offset_source: offset to start copying from
            length: number of entries to copy
        """
        for name, _ in self._fields:
            this_field = getattr(self, name)

            for j in range(length):
                this_field[offset_target + j] = this_field[offset_source + j]

    def sort(self, compare_function: Callable[[int, int], int]):
        """
        Sort entries in the store given the compare function

        Args:
            compare_function: function to sort by (i, j) -> int
                the return value is negative if item[i] is smaller than item[j]
        """
        this_store = self
        tmp_store = Store(self._fields, 1)

        def swap(index1: int, index2: int):
            if index1 == index2:
                return
            tmp_store.copy_from(this_store, 0, index1, 1)
            this_store.copy_within(index1, index2, 1)
            this_store.copy_from(tmp_store, index2, 0, 1)

        def quicksort(left: int, right: int):
            if left < right:
                pivot = (left + right) // 2
                left_new = left
                right_new = right

                while True:
                    while compare_function(left_new, pivot) < 0:
                        left_new += 1

                    while compare_function(right_new, pivot) > 0:
                        right_new -= 1

                    if left_new <= right_new:
                        if left_new == pivot:
                            pivot = right_new
                        elif right_new == pivot:
                            pivot = left_new

                        swap(left_new, right_new)
                        left_new += 1
                        right_new -= 1

                    if left_new > right_new:
                        break

                quicksort(left, right_new)
                quicksort(left_new, right)

        quicksort(0, self.count - 1)

    def clear(self):
        """
        Empty the store, but keep the allocated memory
        """
        self.count = 0

    def dispose(self):
        """
        Dispose of the store entries and fields
        """
        del self.capacity
        del self.count

        for name, _ in self._fields:
            if hasattr(self, name):
                delattr(self, name)
