import math
import os
import random
import sys
import unittest

# Add src directory to Python path
this_dir = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(this_dir, ".."))

from src.pdbstruct import vector3d


def get_random_real():
    """Returns a random number from [-1.0, 1.0]"""
    return (random.random() - 0.5) * 2.0


def get_random_vector():
    """Returns a random vector in the unit circle"""
    return vector3d.Vector3d(get_random_real(), get_random_real(), get_random_real())


def radians(degrees):
    """Convert degrees to radians"""
    return (degrees / 180.0) * math.pi


def random_vec():
    return vector3d.Vector3d(
        random.uniform(-100, 100),
        random.uniform(-100, 100),
        random.uniform(-100, 100),
    )


def random_origin_rotation():
    axis = random_vec()
    angle = random.uniform(-math.pi / 2.0, math.pi / 2.0)
    return vector3d.rotation_at_origin(axis, angle)


def random_transform():
    axis = random_vec()
    angle = random.uniform(-math.pi / 2.0, math.pi / 2.0)
    center = random_vec()
    return vector3d.rotation_at_center(axis, angle, center)


class TestVector3d(unittest.TestCase):
    """Test cases for Vector3d class and related functions."""

    def setUp(self):
        """Set up test fixtures before each test method."""
        self.v1 = vector3d.Vector3d(1.0, 2.0, 3.0)
        self.v2 = vector3d.Vector3d(4.0, 5.0, 6.0)
        self.v3 = vector3d.Vector3d(0.0, 0.0, 0.0)
        self.unit_x = vector3d.Vector3d(1.0, 0.0, 0.0)
        self.unit_y = vector3d.Vector3d(0.0, 1.0, 0.0)
        self.unit_z = vector3d.Vector3d(0.0, 0.0, 1.0)

    def test_vector3d_init(self):
        """Test Vector3d initialization."""
        # Default initialization
        v = vector3d.Vector3d()
        self.assertEqual(v.x, 0.0)
        self.assertEqual(v.y, 0.0)
        self.assertEqual(v.z, 0.0)

        # Parameterized initialization
        v = vector3d.Vector3d(1.0, 2.0, 3.0)
        self.assertEqual(v.x, 1.0)
        self.assertEqual(v.y, 2.0)
        self.assertEqual(v.z, 3.0)

    def test_vector3d_indexing(self):
        """Test Vector3d indexing operations."""
        v = vector3d.Vector3d(1.0, 2.0, 3.0)

        # Test getting values
        self.assertEqual(v[0], 1.0)
        self.assertEqual(v[1], 2.0)
        self.assertEqual(v[2], 3.0)

        # Test setting values
        v[0] = 10.0
        v[1] = 20.0
        v[2] = 30.0
        self.assertEqual(v.x, 10.0)
        self.assertEqual(v.y, 20.0)
        self.assertEqual(v.z, 30.0)

        # Test index out of range
        with self.assertRaises(IndexError):
            _ = v[3]
        with self.assertRaises(IndexError):
            v[3] = 1.0

    def test_vector3d_len_and_iter(self):
        """Test Vector3d length and iteration."""
        v = vector3d.Vector3d(1.0, 2.0, 3.0)

        # Test length
        self.assertEqual(len(v), 3)

        # Test iteration
        components = list(v)
        self.assertEqual(components, [1.0, 2.0, 3.0])

    def test_vector3d_arithmetic_operators(self):
        """Test Vector3d arithmetic operators."""
        v1 = vector3d.Vector3d(1.0, 2.0, 3.0)
        v2 = vector3d.Vector3d(4.0, 5.0, 6.0)

        # Addition
        v3 = v1 + v2
        self.assertEqual(v3.x, 5.0)
        self.assertEqual(v3.y, 7.0)
        self.assertEqual(v3.z, 9.0)

        # Subtraction
        v4 = v2 - v1
        self.assertEqual(v4.x, 3.0)
        self.assertEqual(v4.y, 3.0)
        self.assertEqual(v4.z, 3.0)

        # Negation
        v5 = -v1
        self.assertEqual(v5.x, -1.0)
        self.assertEqual(v5.y, -2.0)
        self.assertEqual(v5.z, -3.0)

        # Positive (copy)
        v6 = +v1
        self.assertEqual(v6.x, 1.0)
        self.assertEqual(v6.y, 2.0)
        self.assertEqual(v6.z, 3.0)
        self.assertIsNot(v6, v1)  # Should be a copy

    def test_vector3d_equality(self):
        """Test Vector3d equality operators."""
        v1 = vector3d.Vector3d(1.0, 2.0, 3.0)
        v2 = vector3d.Vector3d(1.0, 2.0, 3.0)
        v3 = vector3d.Vector3d(1.0, 2.0, 3.1)

        self.assertTrue(v1 == v2)
        self.assertFalse(v1 == v3)
        self.assertFalse(v1 != v2)
        self.assertTrue(v1 != v3)

    def test_vector3d_string_representation(self):
        """Test Vector3d string representations."""
        v = vector3d.Vector3d(1.234, 2.567, 3.891)

        str_repr = str(v)
        self.assertIn("1.23", str_repr)
        self.assertIn("2.57", str_repr)
        self.assertIn("3.89", str_repr)

        repr_str = repr(v)
        self.assertIn("Vector3d", repr_str)
        self.assertIn("1.234", repr_str)

    def test_vector3d_methods(self):
        """Test Vector3d instance methods."""
        v = vector3d.Vector3d(3.0, 4.0, 0.0)

        # Test set method
        v.set(1.0, 2.0, 3.0)
        self.assertEqual(v.x, 1.0)
        self.assertEqual(v.y, 2.0)
        self.assertEqual(v.z, 3.0)

        # Test copy method
        v_copy = v.copy()
        self.assertEqual(v_copy.x, v.x)
        self.assertEqual(v_copy.y, v.y)
        self.assertEqual(v_copy.z, v.z)
        self.assertIsNot(v_copy, v)

        # Test length methods
        v.set(3.0, 4.0, 0.0)
        self.assertEqual(v.mag(), 5.0)

        # Test scale method
        v.scale(2.0)
        self.assertEqual(v.x, 6.0)
        self.assertEqual(v.y, 8.0)
        self.assertEqual(v.z, 0.0)

        # Test scaled_vec method
        v.set(1.0, 2.0, 3.0)
        w = vector3d.scaled_vec(v, 2.0)
        self.assertEqual(w.x, 2.0)
        self.assertEqual(w.y, 4.0)
        self.assertEqual(w.z, 6.0)
        # Original should be unchanged
        self.assertEqual(v.x, 1.0)

        # Test normalize method
        v.set(3.0, 4.0, 0.0)
        v.normalize()
        self.assertAlmostEqual(v.mag(), 1.0, places=6)

        # Test normal_vec method
        v.set(3.0, 4.0, 0.0)
        v.normalize()
        self.assertAlmostEqual(v.mag(), 1.0, places=6)

        # Test tuple method
        v.set(1.0, 2.0, 3.0)
        t = v.tuple()
        self.assertEqual(t, (1.0, 2.0, 3.0))

    def test_vector_creation_and_equality(self):
        """Test vector creation and cloning"""
        v = get_random_vector()
        w = v.copy()
        self.assertTrue(vector3d.is_vec_equal(v, w))

    def test_vector_addition_subtraction(self):
        """Test vector additions and subtractions"""
        v = get_random_vector()
        b = get_random_vector()
        c = vector3d.add_vec(v, b)
        d = vector3d.sub_vec(c, b)
        self.assertTrue(vector3d.is_vec_equal(v, d))

    def test_vector_length_functions(self):
        """Test vector length calculation functions."""
        v = vector3d.Vector3d(3.0, 4.0, 0.0)

        self.assertEqual(vector3d.vec_length_sq(v), 25.0)
        self.assertEqual(vector3d.vec_length(v), 5.0)

    def test_vector_scaling_functions(self):
        """Test vector scaling functions."""
        v = vector3d.Vector3d(1.0, 2.0, 3.0)

        # Test scale_vec (in-place)
        v.scale(2.0)
        self.assertEqual(v.x, 2.0)
        self.assertEqual(v.y, 4.0)
        self.assertEqual(v.z, 6.0)

        # Test scale_vec_copy
        v.set(1.0, 2.0, 3.0)
        v_scaled = vector3d.scaled_vec(v, 3.0)
        self.assertEqual(v_scaled.x, 3.0)
        self.assertEqual(v_scaled.y, 6.0)
        self.assertEqual(v_scaled.z, 9.0)
        # Original unchanged
        self.assertEqual(v.x, 1.0)

    def test_vector_normalization_functions(self):
        """Test vector normalization functions."""
        v = vector3d.Vector3d(3.0, 4.0, 0.0)

        # Test normalize_vec_inplace
        v.normalize()
        self.assertAlmostEqual(v.mag(), 1.0, places=6)

        # Test normalize_vec
        v.set(6.0, 8.0, 0.0)
        v_norm = vector3d.normalized_vec(v)
        self.assertAlmostEqual(vector3d.vec_length(v_norm), 1.0, places=6)
        # Original unchanged
        self.assertEqual(vector3d.vec_length(v), 10.0)

        # Test normalization of zero vector
        zero_vec = vector3d.Vector3d(0.0, 0.0, 0.0)
        norm_zero = vector3d.normalized_vec(zero_vec)
        self.assertEqual(norm_zero.x, 0.0)
        self.assertEqual(norm_zero.y, 0.0)
        self.assertEqual(norm_zero.z, 0.0)

    def test_parallel_perpendicular_functions(self):
        """Test parallel and perpendicular vector functions."""
        v = vector3d.Vector3d(1.0, 1.0, 0.0)
        axis = vector3d.Vector3d(1.0, 0.0, 0.0)

        # Test parallel_vec
        v_parallel = vector3d.parallel_vec(v, axis)
        self.assertEqual(v_parallel.x, 1.0)
        self.assertAlmostEqual(v_parallel.y, 0.0, places=6)
        self.assertAlmostEqual(v_parallel.z, 0.0, places=6)

        # Test perpendicular_vec
        v_perp = vector3d.perpendicular_vec(v, axis)
        self.assertAlmostEqual(v_perp.x, 0.0, places=6)
        self.assertEqual(v_perp.y, 1.0)
        self.assertEqual(v_perp.z, 0.0)

    def test_dot_product(self):
        """Test dot product function."""
        v1 = vector3d.Vector3d(1.0, 2.0, 3.0)
        v2 = vector3d.Vector3d(4.0, 5.0, 6.0)

        result = vector3d.dot(v1, v2)
        expected = 1.0 * 4.0 + 2.0 * 5.0 + 3.0 * 6.0  # 4 + 10 + 18 = 32
        self.assertEqual(result, expected)

    def test_cross_product(self):
        """Test cross product functions."""
        v1 = vector3d.Vector3d(1.0, 0.0, 0.0)
        v2 = vector3d.Vector3d(0.0, 1.0, 0.0)

        # Test cross_product_vec
        cross = vector3d.cross_product_vec(v1, v2)
        self.assertEqual(cross.x, 0.0)
        self.assertEqual(cross.y, 0.0)
        self.assertEqual(cross.z, 1.0)

    def test_distance_functions(self):
        """Test position distance functions."""
        p1 = vector3d.Vector3d(0.0, 0.0, 0.0)
        p2 = vector3d.Vector3d(3.0, 4.0, 0.0)

        # Test pos_distance_sq
        dist_sq = vector3d.pos_distance_sq(p1, p2)
        self.assertEqual(dist_sq, 25.0)

        # Test pos_distance
        dist = vector3d.pos_distance(p1, p2)
        self.assertEqual(dist, 5.0)

    def test_angle_functions(self):
        """Test angle calculation functions."""
        # Test vec_angle
        v1 = vector3d.Vector3d(1.0, 0.0, 0.0)
        v2 = vector3d.Vector3d(0.0, 1.0, 0.0)

        angle = vector3d.vec_angle(v1, v2)
        self.assertAlmostEqual(angle, math.pi / 2, places=6)

        # Test parallel vectors
        v3 = vector3d.Vector3d(2.0, 0.0, 0.0)
        angle_parallel = vector3d.vec_angle(v1, v3)
        self.assertAlmostEqual(angle_parallel, 0.0, places=6)

        # Test opposite vectors
        v4 = vector3d.Vector3d(-1.0, 0.0, 0.0)
        angle_opposite = vector3d.vec_angle(v1, v4)
        self.assertAlmostEqual(angle_opposite, math.pi, places=6)

        # Test pos_angle
        p1 = vector3d.Vector3d(1.0, 0.0, 0.0)
        p2 = vector3d.Vector3d(0.0, 0.0, 0.0)  # vertex
        p3 = vector3d.Vector3d(0.0, 1.0, 0.0)

        pos_angle_result = vector3d.pos_angle(p1, p2, p3)
        self.assertAlmostEqual(pos_angle_result, math.pi / 2, places=6)

    def test_dihedral_functions(self):
        """Test dihedral angle functions."""
        # Test vec_dihedral
        a = vector3d.Vector3d(1.0, 0.0, 0.0)
        axis = vector3d.Vector3d(0.0, 0.0, 1.0)
        c = vector3d.Vector3d(0.0, 1.0, 0.0)

        dihedral = vector3d.vec_dihedral(a, axis, c)
        self.assertAlmostEqual(abs(dihedral), math.pi / 2, places=6)

        # Test pos_dihedral
        p1 = vector3d.Vector3d(1.0, 0.0, 0.0)
        p2 = vector3d.Vector3d(0.0, 0.0, 0.0)
        p3 = vector3d.Vector3d(0.0, 0.0, 1.0)
        p4 = vector3d.Vector3d(0.0, 1.0, 1.0)

        pos_dihedral_result = vector3d.pos_dihedral(p1, p2, p3, p4)
        self.assertAlmostEqual(abs(pos_dihedral_result), math.pi / 2, places=6)

    def test_angle_normalization(self):
        """Test angle normalization functions."""
        # Test normalize_angle
        angle1 = vector3d.normalize_angle(3 * math.pi)
        self.assertAlmostEqual(angle1, math.pi, places=6)

        angle2 = vector3d.normalize_angle(-3 * math.pi)
        self.assertAlmostEqual(angle2, math.pi, places=6)

        angle3 = vector3d.normalize_angle(math.pi / 2)
        self.assertAlmostEqual(angle3, math.pi / 2, places=6)

        # Test angle_diff
        diff = vector3d.angle_diff(math.pi / 4, 3 * math.pi / 4)
        self.assertAlmostEqual(diff, -math.pi / 2, places=6)

    def test_is_near_zero(self):
        """Test is_near_zero utility function."""
        self.assertTrue(vector3d.is_near_zero(1e-7))
        self.assertTrue(vector3d.is_near_zero(0.0))
        self.assertFalse(vector3d.is_near_zero(1e-5))
        self.assertFalse(vector3d.is_near_zero(1.0))

    def test_vec_equal_with_tolerance(self):
        """Test vec_equal function with tolerance."""
        v1 = vector3d.Vector3d(1.0, 2.0, 3.0)
        v2 = vector3d.Vector3d(1.0 + 1e-7, 2.0 + 1e-7, 3.0 + 1e-7)
        v3 = vector3d.Vector3d(1.1, 2.0, 3.0)

        self.assertTrue(vector3d.is_vec_equal(v1, v2))  # Within tolerance
        self.assertFalse(vector3d.is_vec_equal(v1, v3))  # Outside tolerance

    def test_vec_to_tuple(self):
        """Test vec_to_tuple function."""
        v = vector3d.Vector3d(1.5, 2.5, 3.5)
        t = v.tuple()
        self.assertEqual(t, (1.5, 2.5, 3.5))
        self.assertIsInstance(t, tuple)

    def test_constants(self):
        """Test module constants."""
        self.assertAlmostEqual(vector3d.RAD2DEG, 180.0 / math.pi, places=10)
        self.assertAlmostEqual(vector3d.DEG2RAD, math.pi / 180.0, places=10)
        self.assertEqual(vector3d.SMALL, 1e-6)


class TestMatrix3d(unittest.TestCase):
    """Test cases for Matrix3d class."""

    def setUp(self):
        """Set up test fixtures before each test method."""
        self.identity = vector3d.Matrix3d()
        self.v1 = vector3d.Vector3d(1.0, 2.0, 3.0)

    def test_matrix3d_init(self):
        """Test Matrix3d initialization creates identity matrix."""
        m = vector3d.Matrix3d()

        # Check diagonal elements
        self.assertEqual(m.elem00, 1.0)
        self.assertEqual(m.elem11, 1.0)
        self.assertEqual(m.elem22, 1.0)
        self.assertEqual(m.elem33, 1.0)

        # Check some off-diagonal elements
        self.assertEqual(m.elem01, 0.0)
        self.assertEqual(m.elem10, 0.0)
        self.assertEqual(m.elem30, 0.0)

    def test_matrix3d_elem_access(self):
        """Test Matrix3d element access."""
        m = vector3d.Matrix3d()

        # Test getting elements
        self.assertEqual(m.elem(0, 0), 1.0)
        self.assertEqual(m.elem(1, 1), 1.0)
        self.assertEqual(m.elem(0, 1), 0.0)

    def test_matrix3d_string_representation(self):
        """Test Matrix3d string representation."""
        m = vector3d.Matrix3d()
        str_repr = str(m)

        self.assertIn("1.00", str_repr)
        self.assertIn("0.00", str_repr)
        self.assertIn("[", str_repr)
        self.assertIn("]", str_repr)

    def test_transform_vec_functions(self):
        """Test vector transformation functions."""
        v = vector3d.Vector3d(1.0, 2.0, 3.0)
        m = vector3d.Matrix3d()  # Identity matrix

        # Test transform_vec_inplace
        v_copy = v.copy()
        v_copy.transform(m)
        self.assertTrue(
            vector3d.is_vec_equal(v, v_copy)
        )  # Should be unchanged by identity

        # Test transform_vec
        v_transformed = vector3d.transformed_vec(v, m)
        self.assertTrue(
            vector3d.is_vec_equal(v, v_transformed)
        )  # Should be unchanged by identity
        self.assertIsNot(v, v_transformed)  # Should be a copy

    def test_matrix3d_transform_vec_method(self):
        """Test Matrix3d transform_vec method."""
        m = vector3d.Matrix3d()  # Identity matrix
        v = vector3d.Vector3d(1.0, 2.0, 3.0)

        v_transformed = m.transformed_vec(v)
        self.assertTrue(vector3d.is_vec_equal(v, v_transformed))
        self.assertIsNot(v, v_transformed)

    def test_orthogonal_rotations(self):
        """Test orthogonal rotations"""
        x = vector3d.Vector3d(random.random(), 0, 0)
        y = vector3d.Vector3d(0, random.random(), 0)
        z = vector3d.Vector3d(0, 0, random.random())

        # Create rotation matrix around y-axis by 90 degrees
        rotation = vector3d.rotation_at_origin(y, radians(90))
        ry_x = rotation.transformed_vec(x)
        ry_x.scale(-1)

        self.assertTrue(vector3d.is_vec_parallel(ry_x, z))

    def test_cross_product(self):
        """Test cross product"""
        x = vector3d.Vector3d(random.random(), 0, 0)
        y = vector3d.Vector3d(0, random.random(), 0)
        z = vector3d.Vector3d(0, 0, random.random())

        cross_x_y = vector3d.cross_product_vec(x, y)
        cross_x_y.normalize()
        z.normalize()
        self.assertTrue(vector3d.is_vec_equal(cross_x_y, z))

        cross_y_x = vector3d.cross_product_vec(y, x)  # Note: corrected order
        neg_z = z.copy()
        neg_z.scale(-1)
        self.assertTrue(vector3d.is_vec_parallel(cross_y_x, neg_z))

    def test_translation(self):
        """Test translation"""
        x = vector3d.Vector3d(get_random_real(), 0, 0)
        y = vector3d.Vector3d(0, get_random_real(), 0)

        translation = vector3d.translation(y)
        x_and_y = translation.transformed_vec(x)
        x_plus_y = vector3d.add_vec(x, y)

        self.assertTrue(vector3d.is_vec_equal(x_plus_y, x_and_y))

    def test_rotation_preserves_length(self):
        """Test that rotation preserves vector length"""
        x = get_random_vector()
        rotation = vector3d.rotation_at_origin(get_random_vector(), random.random())
        y = rotation.transformed_vec(x)
        self.assertAlmostEqual(x.mag(), y.mag(), 6)

    def test_rotated_pos(self):
        """Test RotatedPos function."""
        # Test rotation of a position around an axis
        theta = math.pi / 2  # 90 degrees
        anchor = vector3d.Vector3d(0.0, 0.0, 0.0)
        center = vector3d.Vector3d(0.0, 0.0, 1.0)
        pos = vector3d.Vector3d(1.0, 0.0, 1.0)

        rotated_pos = vector3d.rotated_pos(theta, anchor, center, pos)
        # This function depends on the Rotation class being properly implemented
        # rotated = RotatedPos(theta, anchor, center, pos)
        # Would need to verify the rotation is correct

    def test_matrix_combination(self):
        """Test matrix combination"""
        # Create some transformation matrices
        rotation1 = vector3d.rotation_at_origin(get_random_vector(), random.random())
        translation1 = vector3d.translation(get_random_vector())
        rotation2 = vector3d.rotation_at_origin(get_random_vector(), random.random())
        translation2 = vector3d.translation(get_random_vector())

        matrices = [rotation1, translation1, rotation2, translation2]

        # Apply matrices individually
        x = get_random_vector()
        result_individual = x.copy()
        for matrix in matrices:
            result_individual = matrix.transformed_vec(result_individual)

        # Combine matrices first, then apply
        combined_matrix = vector3d.Matrix3d()  # Identity matrix
        for matrix in matrices:
            combined_matrix = matrix * combined_matrix

        result_combined = combined_matrix.transformed_vec(x)

        self.assertTrue(vector3d.is_vec_equal(result_individual, result_combined))


class TestAdvancedFunctions(unittest.TestCase):
    """Test cases for advanced geometric functions."""

    def test_edge_cases(self):
        """Test edge cases and error conditions."""
        # Test with very small vectors
        tiny_vec = vector3d.Vector3d(1e-10, 1e-10, 1e-10)
        length = vector3d.vec_length(tiny_vec)
        self.assertGreater(length, 0.0)

        # Test normalization of tiny vector
        normalized = vector3d.normalized_vec(tiny_vec)
        # Should handle gracefully without division by zero

        # Test angle between nearly parallel vectors
        v1 = vector3d.Vector3d(1.0, 0.0, 0.0)
        v2 = vector3d.Vector3d(1.0, 1e-10, 0.0)
        angle = vector3d.vec_angle(v1, v2)
        self.assertAlmostEqual(angle, 0.0, places=6)

    def test_vector_operations_consistency(self):
        """Test consistency between different ways of doing the same operation."""
        v1 = vector3d.Vector3d(1.0, 2.0, 3.0)
        v2 = vector3d.Vector3d(4.0, 5.0, 6.0)

        # Test that operator and function give same result
        sum_op = v1 + v2
        sum_func = vector3d.add_vec(v1, v2)
        self.assertTrue(vector3d.is_vec_equal(sum_op, sum_func))

        diff_op = v1 - v2
        diff_func = vector3d.sub_vec(v1, v2)
        self.assertTrue(vector3d.is_vec_equal(diff_op, diff_func))

        neg_op = -v1
        sum_vec = neg_op + v1
        self.assertTrue(vector3d.is_near_zero(sum_vec.mag()))

    def test_mathematical_properties(self):
        """Test mathematical properties of vector operations."""
        v1 = vector3d.Vector3d(1.0, 2.0, 3.0)
        v2 = vector3d.Vector3d(4.0, 5.0, 6.0)
        v3 = vector3d.Vector3d(7.0, 8.0, 9.0)

        # Test commutativity of addition
        sum1 = vector3d.add_vec(v1, v2)
        sum2 = vector3d.add_vec(v2, v1)
        self.assertTrue(vector3d.is_vec_equal(sum1, sum2))

        # Test associativity of addition
        sum_assoc1 = vector3d.add_vec(vector3d.add_vec(v1, v2), v3)
        sum_assoc2 = vector3d.add_vec(v1, vector3d.add_vec(v2, v3))
        self.assertTrue(vector3d.is_vec_equal(sum_assoc1, sum_assoc2))

        # Test commutativity of dot product
        dot1 = vector3d.dot(v1, v2)
        dot2 = vector3d.dot(v2, v1)
        self.assertEqual(dot1, dot2)

        # Test anti-commutativity of cross product
        cross1 = vector3d.cross_product_vec(v1, v2)
        cross2 = vector3d.cross_product_vec(v2, v1)
        cross2_neg = -cross2
        self.assertTrue(vector3d.is_vec_equal(cross1, cross2_neg))

    def test_normalization_properties(self):
        """Test properties of vector normalization."""
        v = vector3d.Vector3d(3.0, 4.0, 5.0)
        original_length = vector3d.vec_length(v)

        # Test that normalized vector has unit length
        v_norm = vector3d.normalized_vec(v)
        norm_length = vector3d.vec_length(v_norm)
        self.assertAlmostEqual(norm_length, 1.0, places=6)

        # Test that direction is preserved (dot product with original)
        dot_product = vector3d.dot(v, v_norm)
        self.assertAlmostEqual(dot_product, original_length, places=6)

    def test_orthogonality_properties(self):
        """Test orthogonality properties."""
        v = vector3d.Vector3d(1.0, 1.0, 0.0)
        axis = vector3d.Vector3d(1.0, 0.0, 0.0)

        # Get parallel and perpendicular components
        v_parallel = vector3d.parallel_vec(v, axis)
        v_perp = vector3d.perpendicular_vec(v, axis)

        # Test that parallel + perpendicular = original
        v_reconstructed = vector3d.add_vec(v_parallel, v_perp)
        self.assertTrue(vector3d.is_vec_equal(v, v_reconstructed))

        # Test that perpendicular component is orthogonal to axis
        dot_perp_axis = vector3d.dot(v_perp, axis)
        self.assertAlmostEqual(dot_perp_axis, 0.0, places=6)

        # Test that parallel component is parallel to axis (cross product = 0)
        cross_parallel_axis = vector3d.cross_product_vec(v_parallel, axis)
        self.assertAlmostEqual(vector3d.vec_length(cross_parallel_axis), 0.0, places=6)

class TestVector3dIntegration(unittest.TestCase):
    """Integration tests for Vector3d operations working together."""

    def test_orthogonal_basis_construction(self):
        """Test construction of orthogonal basis from arbitrary vector."""
        # Start with arbitrary vector
        v1 = vector3d.Vector3d(1.0, 2.0, 3.0)
        v1_norm = vector3d.normalized_vec(v1)

        # Create second vector not parallel to first
        temp = vector3d.Vector3d(0.0, 1.0, 0.0)
        if abs(vector3d.dot(v1_norm, temp)) > 0.9:  # Nearly parallel
            temp = vector3d.Vector3d(1.0, 0.0, 0.0)

        # Create orthogonal basis
        v2 = vector3d.perpendicular_vec(temp, v1_norm)
        v2_norm = vector3d.normalized_vec(v2)

        v3 = vector3d.cross_product_vec(v1_norm, v2_norm)
        v3_norm = vector3d.normalized_vec(v3)

        # Test orthogonality
        self.assertAlmostEqual(vector3d.dot(v1_norm, v2_norm), 0.0, places=6)
        self.assertAlmostEqual(vector3d.dot(v1_norm, v3_norm), 0.0, places=6)
        self.assertAlmostEqual(vector3d.dot(v2_norm, v3_norm), 0.0, places=6)

        # Test unit length
        self.assertAlmostEqual(vector3d.vec_length(v1_norm), 1.0, places=6)
        self.assertAlmostEqual(vector3d.vec_length(v2_norm), 1.0, places=6)
        self.assertAlmostEqual(vector3d.vec_length(v3_norm), 1.0, places=6)

    def test_triangle_properties(self):
        """Test geometric properties of triangles using vector operations."""
        # Define triangle vertices
        a = vector3d.Vector3d(0.0, 0.0, 0.0)
        b = vector3d.Vector3d(3.0, 0.0, 0.0)
        c = vector3d.Vector3d(0.0, 4.0, 0.0)

        # Calculate side lengths using vectors
        ab = b - a
        bc = c - b
        ca = a - c

        side_a = vector3d.vec_length(bc)  # Opposite to vertex a
        side_b = vector3d.vec_length(ca)  # Opposite to vertex b
        side_c = vector3d.vec_length(ab)  # Opposite to vertex c

        # This is a 3-4-5 right triangle
        self.assertAlmostEqual(side_a, 5.0, places=6)
        self.assertAlmostEqual(side_b, 4.0, places=6)
        self.assertAlmostEqual(side_c, 3.0, places=6)

        # Test angles
        angle_at_a = vector3d.pos_angle(b, a, c)
        angle_at_b = vector3d.pos_angle(a, b, c)
        angle_at_c = vector3d.pos_angle(a, c, b)

        # Sum of angles should be π
        angle_sum = angle_at_a + angle_at_b + angle_at_c
        self.assertAlmostEqual(angle_sum, math.pi, places=6)

        # Right angle at origin
        self.assertAlmostEqual(angle_at_a, math.pi / 2, places=6)

    def test_vector_projection_decomposition(self):
        """Test vector projection and decomposition."""
        # Original vector
        v = vector3d.Vector3d(3.0, 4.0, 5.0)

        # Project onto coordinate axes
        x_axis = vector3d.Vector3d(1.0, 0.0, 0.0)
        y_axis = vector3d.Vector3d(0.0, 1.0, 0.0)
        z_axis = vector3d.Vector3d(0.0, 0.0, 1.0)

        proj_x = vector3d.parallel_vec(v, x_axis)
        proj_y = vector3d.parallel_vec(v, y_axis)
        proj_z = vector3d.parallel_vec(v, z_axis)

        # Projections should sum to original vector
        reconstructed = proj_x + proj_y + proj_z

        self.assertAlmostEqual(reconstructed.x, v.x, places=6)
        self.assertAlmostEqual(reconstructed.y, v.y, places=6)
        self.assertAlmostEqual(reconstructed.z, v.z, places=6)

    def test_dihedral_angle_calculation(self):
        """Test dihedral angle calculation in molecular context."""
        # Define four points representing a dihedral angle
        # Like in a protein backbone: N-CA-CB-CG
        p1 = vector3d.Vector3d(0.0, 0.0, 0.0)  # N
        p2 = vector3d.Vector3d(1.0, 0.0, 0.0)  # CA
        p3 = vector3d.Vector3d(1.0, 1.0, 0.0)  # CB
        p4 = vector3d.Vector3d(2.0, 1.0, 0.0)  # CG (trans)

        dihedral = vector3d.pos_dihedral(p1, p2, p3, p4)

        # This should be a 180° (π radians) dihedral angle
        self.assertAlmostEqual(abs(dihedral), math.pi, places=6)

        # Test with cis configuration
        p4_cis = vector3d.Vector3d(0.0, 1.0, 0.0)
        dihedral_cis = vector3d.pos_dihedral(p1, p2, p3, p4_cis)

        # This should be close to 0° dihedral angle
        self.assertAlmostEqual(abs(dihedral_cis), 0.0, places=6)

    def test_coordinate_system_transformation(self):
        """Test transformation between coordinate systems."""
        # Define a point in one coordinate system
        point = vector3d.Vector3d(1.0, 2.0, 3.0)

        # Identity transformation should leave point unchanged
        identity = vector3d.Matrix3d()
        transformed = vector3d.transformed_vec(point, identity)

        self.assertAlmostEqual(transformed.x, point.x, places=6)
        self.assertAlmostEqual(transformed.y, point.y, places=6)
        self.assertAlmostEqual(transformed.z, point.z, places=6)

if __name__ == "__main__":
    unittest.main()
