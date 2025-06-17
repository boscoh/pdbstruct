import math

RAD2DEG = 180.0 / math.pi
DEG2RAD = math.pi / 180.0
SMALL = 1e-6


def is_near_zero(a):
    return a < SMALL


class Vector3d:
    """
    A 3D vector class for geometric calculations in 3D space.

    This class represents a vector in 3D Cartesian coordinates and provides
    methods for common vector operations including arithmetic, transformations,
    and geometric calculations.

    Attributes:
        x (float): The x-coordinate component of the vector
        y (float): The y-coordinate component of the vector
        z (float): The z-coordinate component of the vector

    Examples:
        >>> v1 = Vector3d(1.0, 2.0, 3.0)
        >>> v2 = Vector3d(4.0, 5.0, 6.0)
        >>> v3 = add_vec(v1, v2)
        >>> print(v3)
        ( 5.00, 7.00, 9.00 )

        >>> length = vec_length(v1)
        >>> normalized = normalized_vec(v1)
        >>> dot_product = dot(v1, v2)

    Note:
        This class supports indexing (v[0], v[1], v[2]), iteration, and
        standard arithmetic operations (+, -, unary -, ==, !=).
        Vectors can be transformed using Matrix3d transformation matrices.
    """

    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def __getitem__(self, index):
        """Allow indexing: v[0] = x, v[1] = y, v[2] = z"""
        if index == 0:
            return self.x
        elif index == 1:
            return self.y
        elif index == 2:
            return self.z
        else:
            raise IndexError("Vector3d index out of range (0-2)")

    def __setitem__(self, index, value):
        """Allow setting via indexing: v[0] = x, v[1] = y, v[2] = z"""
        if index == 0:
            self.x = value
        elif index == 1:
            self.y = value
        elif index == 2:
            self.z = value
        else:
            raise IndexError("Vector3d index out of range (0-2)")

    def __len__(self):
        """Return length of vector (always 3 for 3D vector)"""
        return 3

    def __iter__(self):
        """Allow iteration over vector components"""
        yield self.x
        yield self.y
        yield self.z

    def __add__(self, rhs):
        return add_vec(self, rhs)

    def __sub__(self, rhs):
        return sub_vec(self, rhs)

    def __neg__(self):
        return Vector3d(-self[0], -self[1], -self[2])

    def __pos__(self):
        return self.copy()

    def __eq__(self, rhs):
        return is_vec_equal(self, rhs)

    def __ne__(self, rhs):
        return not is_vec_equal(self, rhs)

    def __str__(self):
        return "( %.2f, %.2f, %.2f )" % (self[0], self[1], self[2])

    def __repr__(self):
        return "Vector3d( %f, %f, %f )" % (self[0], self[1], self[2])

    def set(self, x, y, z):
        """Set vector components"""
        self[0] = x
        self[1] = y
        self[2] = z

    def copy(self):
        return Vector3d(self[0], self[1], self[2])

    def mag(self):
        return vec_length(self)

    def scale(self, scale):
        self[0] *= scale
        self[1] *= scale
        self[2] *= scale

    def normalize(self):
        length = vec_length(self)
        if not is_near_zero(length):
            self.scale(1.0 / length)

    def transform(self, matrix):
        x = (
            matrix.elem00 * self[0]
            + matrix.elem10 * self[1]
            + matrix.elem20 * self[2]
            + matrix.elem30
        )
        y = (
            matrix.elem01 * self[0]
            + matrix.elem11 * self[1]
            + matrix.elem21 * self[2]
            + matrix.elem31
        )
        z = (
            matrix.elem02 * self[0]
            + matrix.elem12 * self[1]
            + matrix.elem22 * self[2]
            + matrix.elem32
        )
        self[0], self[1], self[2] = x, y, z

    def tuple(self):
        return (self[0], self[1], self[2])


def add_vec(a, b):
    """Add two vectors"""
    return Vector3d(a[0] + b[0], a[1] + b[1], a[2] + b[2])


def sub_vec(a, b):
    """Subtract vector b from vector a"""
    return Vector3d(a[0] - b[0], a[1] - b[1], a[2] - b[2])


def is_vec_equal(a, b):
    """Check if two vectors are equal within tolerance"""
    return (
        is_near_zero(math.fabs(a[0] - b[0]))
        and is_near_zero(math.fabs(a[1] - b[1]))
        and is_near_zero(math.fabs(a[2] - b[2]))
    )


def vec_length_sq(v):
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2]


def vec_length(v):
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def pos_distance(p1, p2):
    return math.sqrt(pos_distance_sq(p1, p2))


def pos_distance_sq(p1, p2):
    x = p1[0] - p2[0]
    y = p1[1] - p2[1]
    z = p1[2] - p2[2]
    return x * x + y * y + z * z


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross_product_vec(a, b):
    """Calculate cross product of two vectors"""
    return Vector3d(
        a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]
    )


def is_vec_parallel(a, b, tolerance=SMALL) -> bool:
    """
    Test if two vectors are parallel (or anti-parallel).
    """
    if vec_length(a) < tolerance or vec_length(b) < tolerance:
        return True
    cross = cross_product_vec(a, b)
    return vec_length(cross) < tolerance


def scaled_vec(v, scale):
    """Return scaled copy of vector"""
    result = v.copy()
    result.scale(scale)
    return result


def normalized_vec(v):
    """Return normalized copy of vector"""
    length = vec_length(v)
    if is_near_zero(length):
        return Vector3d(0.0, 0.0, 0.0)
    return scaled_vec(v, 1.0 / length)


def parallel_vec(v, axis):
    """Return component of vector parallel to axis"""
    axis_len = vec_length(axis)
    if is_near_zero(axis_len):
        return v.copy()
    else:
        dot_product = dot(v, axis)
        axis_len_sq = axis_len * axis_len
        return scaled_vec(axis, dot_product / axis_len_sq)


def perpendicular_vec(v, axis):
    """Return component of vector perpendicular to axis"""
    parallel = parallel_vec(v, axis)
    return sub_vec(v, parallel)


def transformed_vec(v, matrix):
    """Return transformed copy of vector"""
    result = v.copy()
    result.transform(matrix)
    return result


def normalize_angle(angle):
    while abs(angle) > math.pi:
        if angle > math.pi:
            angle -= math.pi * 2
        if angle < -math.pi:
            angle += 2 * math.pi
    if is_near_zero(abs(angle + math.pi)):
        angle = math.pi
    return angle


def angle_diff(angle1, angle2):
    norm_angle1 = normalize_angle(angle1)
    norm_angle2 = normalize_angle(angle2)
    return normalize_angle(norm_angle1 - norm_angle2)


def vec_angle(a, b):
    a_len = vec_length(a)
    b_len = vec_length(b)

    if a_len * b_len < 1e-6:
        return 0.0

    c = 1.0 * dot(a, b) / a_len / b_len

    if c >= 1.0:
        return 0.0
    elif c <= -1.0:
        return math.pi
    else:
        return math.acos(c)


def pos_angle(p1, p2, p3):
    return vec_angle(sub_vec(p1, p2), sub_vec(p3, p2))


def vec_dihedral(a, axis, c):
    ap = perpendicular_vec(a, axis)
    cp = perpendicular_vec(c, axis)

    angle = vec_angle(ap, cp)

    if dot(cross_product_vec(ap, cp), axis) > 0:
        angle = -angle

    return angle


def pos_dihedral(p1, p2, p3, p4):
    return vec_dihedral(sub_vec(p1, p2), sub_vec(p2, p3), sub_vec(p4, p3))


def rotated_pos(theta, anchor, center, pos):
    return rotation_at_center(sub_vec(center, anchor), theta, center).transformed_vec(pos)


class Matrix3d:
    def __init__(self):
        self.elem00 = 1.0
        self.elem01 = 0.0
        self.elem02 = 0.0
        self.elem03 = 0.0

        self.elem10 = 0.0
        self.elem11 = 1.0
        self.elem12 = 0.0
        self.elem13 = 0.0

        self.elem20 = 0.0
        self.elem21 = 0.0
        self.elem22 = 1.0
        self.elem23 = 0.0

        self.elem30 = 0.0
        self.elem31 = 0.0
        self.elem32 = 0.0
        self.elem33 = 1.0

    def __str__(self):
        row1 = "  [% .2f, % .2f, % .2f ]\n" % (self.elem00, self.elem01, self.elem02)
        row2 = "  [% .2f, % .2f, % .2f ]\n" % (self.elem10, self.elem11, self.elem12)
        row3 = "  [% .2f, % .2f, % .2f ]\n" % (self.elem20, self.elem21, self.elem22)
        row4 = "  [ ------------------ ]\n"
        row5 = "  [% .2f, % .2f, % .2f ]" % (self.elem30, self.elem31, self.elem32)
        return row1 + row2 + row3 + row4 + row5

    def elem(self, i, j):
        if j == 0:
            if i == 0:
                return self.elem00
            if i == 1:
                return self.elem10
            if i == 2:
                return self.elem20
            if i == 3:
                return self.elem30
        if j == 1:
            if i == 0:
                return self.elem01
            if i == 1:
                return self.elem11
            if i == 2:
                return self.elem21
            if i == 3:
                return self.elem31
        if j == 2:
            if i == 0:
                return self.elem02
            if i == 1:
                return self.elem12
            if i == 2:
                return self.elem22
            if i == 3:
                return self.elem32
        if j == 3:
            if i == 0:
                return self.elem03
            if i == 1:
                return self.elem13
            if i == 2:
                return self.elem23
            if i == 3:
                return self.elem33

    def set_elem(self, i, j, val):
        if j == 0:
            if i == 0:
                self.elem00 = val
            if i == 1:
                self.elem10 = val
            if i == 2:
                self.elem20 = val
            if i == 3:
                self.elem30 = val
        if j == 1:
            if i == 0:
                self.elem01 = val
            if i == 1:
                self.elem11 = val
            if i == 2:
                self.elem21 = val
            if i == 3:
                self.elem31 = val
        if j == 2:
            if i == 0:
                self.elem02 = val
            if i == 1:
                self.elem12 = val
            if i == 2:
                self.elem22 = val
            if i == 3:
                self.elem32 = val
        if j == 3:
            if i == 0:
                self.elem03 = val
            if i == 1:
                self.elem13 = val
            if i == 2:
                self.elem23 = val
            if i == 3:
                self.elem33 = val

    def __eq__(self, rhs):
        for i in range(0, 3):
            for j in range(0, 3):
                if abs(self.elem(i, j) - rhs.elem(i, j)) > SMALL:
                    return False
        return True

    def __mul__(self, rhs):
        c = Matrix3d()
        for i in range(0, 3):
            for j in range(0, 3):
                val = 0.0
                for k in range(0, 3):
                    val += self.elem(k, i) * rhs.elem(j, k)
                c.set_elem(j, i, val)
            # c(3,i) is the translation vector
            val = self.elem(3, i)
            for k in range(0, 3):
                val += self.elem(k, i) * rhs.elem(3, k)
            c.set_elem(3, i, val)
        return c

    def transformed_vec(self, v):
        # v'[i] = sum(over j) M[j][i] v[j]
        x = self.elem00 * v[0] + self.elem10 * v[1] + self.elem20 * v[2] + self.elem30
        y = self.elem01 * v[0] + self.elem11 * v[1] + self.elem21 * v[2] + self.elem31
        z = self.elem02 * v[0] + self.elem12 * v[1] + self.elem22 * v[2] + self.elem32
        return Vector3d(x, y, z)


def rotation_at_origin(axis, theta):
    """matrix to rotate a vector at origin"""
    v = normalized_vec(axis)

    c = math.cos(float(theta))
    s = math.sin(float(theta))
    t = 1.0 - c

    m = Matrix3d()

    m.elem00 = t * v[0] * v[0] + c
    m.elem01 = t * v[0] * v[1] + v[2] * s
    m.elem02 = t * v[0] * v[2] - v[1] * s

    m.elem10 = t * v[1] * v[0] - v[2] * s
    m.elem11 = t * v[1] * v[1] + c
    m.elem12 = t * v[1] * v[2] + v[0] * s

    m.elem20 = t * v[2] * v[0] + v[1] * s
    m.elem21 = t * v[2] * v[1] - v[0] * s
    m.elem22 = t * v[2] * v[2] + c

    return m


def translation(p):
    """matrix to translate a vector"""
    m = Matrix3d()
    m.elem30 = p[0]
    m.elem31 = p[1]
    m.elem32 = p[2]
    return m


def rotation_at_center(axis, theta, center):
    """matrix to rotate around an axis at center"""
    rot = rotation_at_origin(axis, theta)
    trans = translation(sub_vec(center, rot.transformed_vec(center)))
    return trans * rot


def superposition(ref1, ref2, ref3, mov1, mov2, mov3):
    mov_diff = sub_vec(mov2, mov1)
    ref_diff = sub_vec(ref2, ref1)

    m1: Matrix3d
    if math.fabs(vec_angle(mov_diff, ref_diff)) < SMALL:
        m1 = translation(sub_vec(ref1, mov1))
    else:
        axis = cross_product_vec(mov_diff, ref_diff)
        torsion = vec_dihedral(ref_diff, axis, mov_diff)
        rot = rotation_at_origin(axis, torsion)
        trans = translation(sub_vec(ref2, rot.transformed_vec(mov2)))
        m1 = trans * rot

    mov_diff = sub_vec(ref2, m1.transformed_vec(mov3))
    ref_diff = sub_vec(ref2, ref3)

    m: Matrix3d
    if math.fabs(vec_angle(mov_diff, ref_diff)) < SMALL:
        m = m1
    else:
        axis = sub_vec(ref2, ref1)
        torsion = vec_dihedral(ref_diff, axis, mov_diff)
        m2 = rotation_at_origin(axis, torsion)
        m3 = translation(sub_vec(ref2, m2.transformed_vec(ref2)))
        m = m3 * m2 * m1

    return m


