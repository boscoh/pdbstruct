import array
import unittest
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


class TestStore(unittest.TestCase):
    """Test cases for the Store class"""

    def setUp(self):
        """Set up test fixtures before each test method."""
        self.fields = [
            ("x", "float32"),  # x coordinate
            ("y", "float32"),  # y coordinate
            ("z", "float32"),  # z coordinate
            ("r", "uint8"),  # red value
            ("g", "uint8"),  # green value
            ("b", "uint8"),  # blue value
            ("a", "uint8"),  # alpha value
            ("id", "int32"),  # id value
        ]

    def test_initialization_with_size(self):
        """Test Store initialization with a specific size"""
        store = Store(self.fields, 10)

        self.assertEqual(store.capacity, 10)
        self.assertEqual(store.count, 0)
        self.assertEqual(len(store._fields), 8)

        # Check if fields are properly initialized
        self.assertTrue(hasattr(store, "x"))
        self.assertTrue(hasattr(store, "y"))
        self.assertTrue(hasattr(store, "z"))
        self.assertTrue(hasattr(store, "r"))
        self.assertTrue(hasattr(store, "g"))
        self.assertTrue(hasattr(store, "b"))
        self.assertTrue(hasattr(store, "a"))
        self.assertTrue(hasattr(store, "id"))

        # Check array sizes - each field has one value per entry
        self.assertEqual(len(store.x), 10)
        self.assertEqual(len(store.y), 10)
        self.assertEqual(len(store.z), 10)
        self.assertEqual(len(store.r), 10)
        self.assertEqual(len(store.g), 10)
        self.assertEqual(len(store.b), 10)
        self.assertEqual(len(store.a), 10)
        self.assertEqual(len(store.id), 10)

    def test_initialization_without_size(self):
        """Test Store initialization without size"""
        store = Store(self.fields)

        self.assertEqual(store.capacity, 0)
        self.assertEqual(store.count, 0)
        self.assertEqual(len(store.x), 0)
        self.assertEqual(len(store.r), 0)
        self.assertEqual(len(store.id), 0)

    def test_add_field(self):
        """Test adding a field to the store"""
        store = Store(self.fields, 5)
        initial_field_count = len(store._fields)

        store.add_field("velocity", "float32")

        self.assertEqual(len(store._fields), initial_field_count + 1)
        self.assertTrue(hasattr(store, "velocity"))
        self.assertEqual(len(store.velocity), 5)

    def test_resize_larger(self):
        """Test resizing store to a larger size"""
        store = Store(self.fields, 5)

        # Set some data
        store.x[0] = 1.0
        store.r[0] = 255
        store.id[0] = 42

        store.resize(10)

        self.assertEqual(store.capacity, 10)
        self.assertEqual(len(store.x), 10)
        self.assertEqual(len(store.r), 10)
        self.assertEqual(len(store.id), 10)

        # Check that data is preserved
        self.assertEqual(store.x[0], 1.0)
        self.assertEqual(store.r[0], 255)
        self.assertEqual(store.id[0], 42)

    def test_resize_smaller(self):
        """Test resizing store to a smaller size"""
        store = Store(self.fields, 10)
        store.count = 8

        # Set some data
        store.x[0] = 1.0
        store.r[0] = 255

        store.resize(5)

        self.assertEqual(store.capacity, 5)
        self.assertEqual(store.count, 5)  # Should be clamped to new length
        self.assertEqual(len(store.x), 5)

        # Check that data is preserved
        self.assertEqual(store.x[0], 1.0)
        self.assertEqual(store.r[0], 255)

    def test_grow_if_full(self):
        """Test automatic growth when store is full"""
        store = Store(self.fields, 2)
        store.count = 2  # Make it full

        store.grow_if_full()

        self.assertEqual(store.capacity, 256)  # 256 if size is too small

    def test_grow_if_full_empty_store(self):
        """Test growth of empty store"""
        store = Store(self.fields, 0)
        store.count = 0

        store.grow_if_full()

        self.assertEqual(store.capacity, 256)  # Should grow to minimum 256

    def test_increment(self):
        """Test increment method"""
        store = Store(self.fields, 2)

        store.increment()
        self.assertEqual(store.count, 1)
        self.assertEqual(store.capacity, 2)

        store.increment()
        self.assertEqual(store.count, 2)
        self.assertEqual(store.capacity, 256)

        # This should trigger growth
        store.increment()
        self.assertEqual(store.count, 3)
        self.assertEqual(store.capacity, 256)

    def test_copy_from(self):
        """Test copying data from another store"""
        store1 = Store(self.fields, 5)
        store2 = Store(self.fields, 5)

        # Set data in store1
        store1.x[0] = 1.0
        store1.y[0] = 2.0
        store1.z[0] = 3.0
        store1.r[0] = 255
        store1.g[0] = 128
        store1.b[0] = 64
        store1.a[0] = 32
        store1.id[0] = 42

        # Copy from store1 to store2
        store2.copy_from(store1, 0, 0, 1)

        # Check that data was copied
        self.assertEqual(store2.x[0], 1.0)
        self.assertEqual(store2.y[0], 2.0)
        self.assertEqual(store2.z[0], 3.0)
        self.assertEqual(store2.r[0], 255)
        self.assertEqual(store2.g[0], 128)
        self.assertEqual(store2.b[0], 64)
        self.assertEqual(store2.a[0], 32)
        self.assertEqual(store2.id[0], 42)

    def test_copy_within(self):
        """Test copying data within the same store"""
        store = Store(self.fields, 5)

        # Set data at index 0
        store.x[0] = 1.0
        store.y[0] = 2.0
        store.z[0] = 3.0
        store.r[0] = 255
        store.id[0] = 42

        # Copy from index 0 to index 1
        store.copy_within(1, 0, 1)

        # Check that data was copied to index 1
        self.assertEqual(store.x[1], 1.0)
        self.assertEqual(store.y[1], 2.0)
        self.assertEqual(store.z[1], 3.0)
        self.assertEqual(store.r[1], 255)
        self.assertEqual(store.id[1], 42)

    def test_sort(self):
        """Test sorting functionality"""
        store = Store([("value", "int32")], 5)
        store.count = 5

        # Set values: [5, 2, 8, 1, 9]
        store.value[0] = 5
        store.value[1] = 2
        store.value[2] = 8
        store.value[3] = 1
        store.value[4] = 9

        # Sort by value (ascending)
        def compare_func(i, j):
            return store.value[i] - store.value[j]

        store.sort(compare_func)

        # Check if sorted: [1, 2, 5, 8, 9]
        expected_values = [1, 2, 5, 8, 9]
        for i in range(5):
            self.assertEqual(store.value[i], expected_values[i])


if __name__ == "__main__":
    # Run tests instead of main when script is executed
    unittest.main()
