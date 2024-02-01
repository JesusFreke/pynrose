#!/usr/bin/python3
# Copyright (c) 2024, Ben Gruver
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
import itertools
import math
import random
from enum import Enum
from functools import cached_property
from typing import List, Optional, Tuple, Sequence


class Vector(object):
    """A 2d vector."""

    def __init__(self, x: float, y: float):
        """Constructs a new (x, y) vector.

        :param x: The component of the vector in the x-axis.
        :param y: The component of the vector in the y-axis.
        """
        self._x = x
        self._y = y

    @property
    def x(self):
        """The component of this vector in the x-axis."""
        return self._x

    @property
    def y(self):
        """The component of this vector in the y-axis."""
        return self._y

    @cached_property
    def length(self):
        """The length of this vector."""
        return math.sqrt(self._x ** 2 + self._y ** 2)

    def dot(self, other: 'Vector'):
        """The dot product between this vector and another vector."""
        return (self._x * other._x) + (self._y * other._y)

    def __eq__(self, other):
        if not isinstance(other, Vector):
            return False
        return self.x == other.x and self.y == other.y

    def __hash__(self):
        return hash((self.x, self.y))

    def __add__(self, other):
        if not isinstance(other, Vector):
            raise TypeError("expecting Vector, got %s" % type(other))
        return Vector(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        if not isinstance(other, Vector):
            raise TypeError("expecting Vector, got %s" % type(other))
        return Vector(self.x - other.x, self.y - other.y)

    def __mul__(self, other):
        if not isinstance(other, int) and not isinstance(other, float):
            raise TypeError("expecting int or float, got %s" % type(other))
        return Vector(self.x * other, self.y * other)

    def __truediv__(self, other):
        if not isinstance(other, int) and not isinstance(other, float):
            raise TypeError("expecting int or float, got %s" % type(other))
        return Vector(self._x / other, self._y / other)

    def __repr__(self):
        return "Vector(%f, %f)" % (self._x, self._y)


class Grid(object):
    """A repeating grid of equal-size grid cells."""

    def __init__(self, origin: Vector, grid_size: Vector):
        """Constructs a new Grid.

        :param origin: The coordinate of the lower left corner of the (0, 0)the grid cell.
        :param grid_size: The size of each grid cell.
        """
        self._origin = origin
        self._grid_size = grid_size

    @property
    def origin(self):
        """The coordinate of the lower left corner of the (0, 0)th grid cell."""
        return self._origin

    @property
    def grid_size(self):
        """The size of each grid cell."""
        return self._grid_size

    def cell(self, x: int, y: int):
        """Returns the (x, y)th grid cell."""
        return GridCell(self, x, y)

    def __eq__(self, other):
        if not isinstance(other, Grid):
            return False
        return self.origin == other.origin and self.grid_size == other.grid_size

    def __hash__(self):
        return hash((self.origin, self.grid_size))


class GridCell(object):
    """An individual grid cell."""

    def __init__(self, grid: Grid, x_multiple: int, y_multiple: int):
        """Constructs a new grid cell at index (x, y) from the given grid definition.

        :param grid: The grid containing this grid cell.
        :param x_multiple: The x multiple of this grid cell.
        :param y_multiple: The y multiple of this grid cell.
        """
        self._grid = grid
        self._x_multiple = x_multiple
        self._y_multiple = y_multiple

        self._origin = grid.origin + Vector(grid.grid_size.x * x_multiple, grid.grid_size.y * y_multiple)
        self._extent = grid.origin + grid.grid_size

    @property
    def grid(self):
        """The grid that this grid cell is a part of."""
        return self._grid

    @property
    def x_multiple(self):
        """The x multiple of this grid cell."""
        return self._x_multiple

    @property
    def y_multiple(self):
        """The y multiple of this grid cell."""
        return self._y_multiple

    @property
    def origin(self):
        """The coordinate of the lower left corner of this grid cell."""
        return self._origin

    @property
    def extent(self):
        """The coordinate of the upper right corner of this grid cell."""
        return self._extent

    def midpoint(self):
        """Returns the midpoint of the grid."""
        return (self.origin + self.extent) / 2

    def corners(self, margin=0.0):
        """Returns the 4 coordinates of the corners of the grid, with an optional margin added.

        :param margin: How much margin to add around the edges of the grid cell.
        """
        yield self.origin + Vector(-margin, -margin)
        yield Vector(self.origin.x - margin, self.extent.y + margin)
        yield self.extent + Vector(margin, margin)
        yield Vector(self.extent.x + margin, self.origin.y - margin)


class RhombusType(Enum):
    THIN = 0
    THICK = 1


class PentAngle(object):
    """A multiple of 360/5 degrees."""

    _SIN = [
        math.sin(0),
        math.sin(1 * 2 * math.pi / 5),
        math.sin(2 * 2 * math.pi / 5),
        math.sin(3 * 2 * math.pi / 5),
        math.sin(4 * 2 * math.pi / 5)]

    _INVERSE_SIN = [
        0,
        1/math.sin(1 * 2 * math.pi / 5),
        1/math.sin(2 * 2 * math.pi / 5),
        1/math.sin(3 * 2 * math.pi / 5),
        1/math.sin(4 * 2 * math.pi / 5)]

    _COS = [
        math.cos(0),
        math.cos(1 * 2 * math.pi / 5),
        math.cos(2 * 2 * math.pi / 5),
        math.cos(3 * 2 * math.pi / 5),
        math.cos(4 * 2 * math.pi / 5)]

    _INVERSE_COS = [
        0,
        1/math.cos(1 * 2 * math.pi / 5),
        1/math.cos(2 * 2 * math.pi / 5),
        1/math.cos(3 * 2 * math.pi / 5),
        1/math.cos(4 * 2 * math.pi / 5)]

    def __init__(self, pentangle: int):
        """Constructs a new PentAngle.

        Prefer using PentAngles.pentangle(pentangle).

        :param pentangle: An integer representing a multiple of 360/5 degrees. Will be normalized to [0, 4].
        """
        self._pentangle: int = int(pentangle % 5)

        self._sin = PentAngle._SIN[pentangle]
        self._cos = PentAngle._COS[pentangle]

    @property
    def pentangle(self):
        """The multiple of 360/5 degrees that this pentangle represents."""
        return self._pentangle

    def degrees(self):
        """The angle of this pentangle in degrees."""
        return self._pentangle * 360 / 5

    def radians(self):
        """The angle of this pentangle in radians."""
        return self._pentangle * 2 * math.pi / 5

    def unit(self) -> 'Vector':
        """Returns a unit vector in the direction of this PentAngle."""
        return Vector(self.sin(), self.cos())

    def sin(self, other: 'PentAngle' = None) -> float:
        """Returns the sin of this pentangle, or of the angle between this pentangle and another pentangle."""
        if other is None:
            return self._sin
        return PentAngle._SIN[(other.pentangle - self.pentangle) % 5]

    def inverse_sin(self, other: 'PentAngle') -> float:
        """Returns the inverse of the sin of the angle between this pentangle and another pentangle."""
        return PentAngle._INVERSE_SIN[(other.pentangle - self.pentangle) % 5]

    def cos(self, other: 'PentAngle' = None) -> float:
        """Returns the cos of this pentangle, or of the angle between this pentangle and another pentangle."""
        if not other:
            return self._cos
        return PentAngle._COS[(other.pentangle - self.pentangle) % 5]

    def inverse_cos(self, other: 'PentAngle') -> float:
        """Returns the inverse of the cos of the angle between this pentangle and another pentangle."""
        return PentAngle._INVERSE_COS[(other.pentangle - self.pentangle) % 5]

    def __eq__(self, other):
        if not isinstance(other, PentAngle):
            return False
        return other.pentangle == self.pentangle

    def __hash__(self):
        return self.pentangle


class PentAngles(object):
    """A convenience class for """

    _ALL = (
        PentAngle(0),
        PentAngle(1),
        PentAngle(2),
        PentAngle(3),
        PentAngle(4)
    )

    @staticmethod
    def pentangle(pentangle: int):
        """Returns a PentAngle object for the given pentangle."""
        return PentAngles._ALL[pentangle]

    @staticmethod
    def all():
        """Returns a tuple of all PentAngles in order."""
        return PentAngles._ALL

    @staticmethod
    def others(pent_angle: 'PentAngle'):
        """Returns the other 4 PentAngles, excluding the given one, in order."""
        for other in PentAngles._ALL:
            if other == pent_angle:
                continue
            yield other


class StripFamily(object):
    """Represents 1 of the 5 families of strips of rhombii in de Bruijn's method."""

    def __init__(self, tiling: 'Tiling', offset: float, pentangle: PentAngle):
        """Constructs a new StripFamily.

        :param tiling: The Tiling that this StripFamily is associated with.
        :param offset: The offset in pentagrid space for this family.
        :param pentangle: Which pentangle this family is associated with.
        """
        self._tiling = tiling
        self._offset = offset
        self._pentangle = pentangle

    @property
    def tiling(self):
        """The Tiling that this StripFamily is associated with."""
        return self._tiling

    @property
    def offset(self):
        """The offset in pentagrid space for this family."""
        return self._offset

    @property
    def pentangle(self):
        """The pentangle for this family."""
        return self._pentangle

    def direction(self) -> Vector:
        """Returns a unit vector in the direction the strips from this family run"""
        return self.pentangle.unit()

    def offset_direction(self) -> Vector:
        """Returns a unit vector perpendicular to the strips from this family, in the positive direction.

        For example, if multiple 0 is a strip that runs vertically at x=0, and multiple 1 is a vertical strip at x=1,
        then offset_direction will be (1, 0).
        """
        direction = self.direction()
        # rotate the vector 90 degrees
        return Vector(direction.y, -direction.x)

    def strip(self, multiple: int):
        """Returns the strip at the given multiple."""
        return Strip(self, multiple)

    def origin(self):
        """Returns the origin of this family."""
        return self.strip(0).origin()

    def strips_near_point(self, point: Vector):
        """Returns 1 or 2 strips from this family that are near the given point.

        The returned strips are ones that could potentially (but may not) contain the rhombus at the given point.
        """
        pentagrid_point = point / 2.5
        multiple = (pentagrid_point - self.origin()).dot(self.offset_direction())

        # the rhombii associated with a strip may wander away from the strip up to a max of just under 2 units in
        # rhombus space, which is .8 units in pentagrid space
        if abs(math.ceil(multiple) - multiple) <= .8:
            yield self.strip(math.ceil(multiple))
        if abs(math.floor(multiple) - multiple) <= .8:
            yield self.strip(math.floor(multiple))


def _det(x1, y1, x2, y2):
    """Returns the determinate of (x1, y1) and (x2, y2)

    | x1 y1 |
    | x2 y2 |
    """
    return x1 * y2 - y1 * x2


def _ccw(p1, p2, p3):
    value = (p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y)
    if math.isclose(value, 0):
        return 0
    return math.copysign(1, value)


def _intersection(p0, p1, p2, p3):
    # see, e.g. https://mathworld.wolfram.com/Line-LineIntersection.html

    denom = _det(
        p0.x - p1.x,
        p0.y - p1.y,
        p2.x - p3.x,
        p2.y - p3.y)

    x_num = _det(
        _det(p0.x, p0.y, p1.x, p1.y),
        p0.x - p1.x,
        _det(p2.x, p2.y, p3.x, p3.y),
        p2.x - p3.x)

    y_num = _det(
        _det(p0.x, p0.y, p1.x, p1.y),
        p0.y - p1.y,
        _det(p2.x, p2.y, p3.x, p3.y),
        p2.y - p3.y)

    return Vector(x_num/denom, y_num/denom)


class Strip(object):
    """Represents a single string of rhombii in de Bruijn's method"""
    def __init__(self, family: StripFamily, multiple: int):
        """Constructs a new Strip.

        :param family: The StringFamily this Strip is associated with.
        :param multiple: The multiple of the Strip.
        """
        self._family = family
        self._multiple = multiple

    @property
    def family(self):
        """The StripFamily this Strip is associated with."""
        return self._family

    @property
    def multiple(self):
        """The multiple of this Strip."""
        return self._multiple

    def rhombii(self, start_distance: float, forward: bool):
        """Infinitely iterates over the rhombii in this Strip.

        :param start_distance: The distance in pentagrid space from this Strip's origin to start iterating.
        :param forward: The direction to iterate.
        """
        return self._rhombii(start_distance, forward)

    def rhombus_at_intersection(self, other: 'Strip'):
        """Returns the Rhombus represented by the intersection of this Strip and the given Strip.

        :param other: The Strip to find the intersection with. Must not be from the same StripFamily as this Strip
                      (i.e. not parallel with this Strip).
        """
        if other.family == self.family:
            raise ValueError("The strips don't intersect")

        lattice_coords = [0] * 5
        lattice_coords[self.family.pentangle.pentangle] = self.multiple

        distance = self.intersection_distance_from_point(other)

        for pentangle in PentAngles.others(self.family.pentangle):
            if pentangle == other.family.pentangle:
                lattice_coords[other.family.pentangle.pentangle] = other.multiple
                continue

            other_family = self.family.tiling.strip_family(pentangle)

            initial_intersection = self.intersection_distance_from_point(other_family.strip(0))
            delta = initial_intersection - distance
            inverse_sin = self.family.pentangle.inverse_sin(pentangle)

            if inverse_sin > 0:
                multiple = int(math.floor(delta / inverse_sin))
                lattice_coords[pentangle.pentangle] = multiple
            else:
                multiple = int(math.ceil(delta / inverse_sin))
                lattice_coords[pentangle.pentangle] = multiple - 1

        return Rhombus(self, other, tuple(lattice_coords))

    def rhombus(self, target_distance: float):
        """Returns the rhombus nearest the given distance from this Strip's origin.

        This finds the intersection closest to the point that is the given distance from this Strip's origin, and
        returns the associated rhombus.

        :param target_distance: The distance from the Strip's origin in pentagrid space.
        """
        closest_value = math.inf
        closest = None

        forward_rhombus = next(self._rhombii(target_distance, True))
        backward_rhombus = next(self._rhombii(target_distance, False))

        for rhombus in (forward_rhombus, backward_rhombus):
            distance_from_target = abs(self.intersection_distance_from_point(rhombus.strip2) - target_distance)
            if distance_from_target < closest_value:
                closest_value = distance_from_target
                closest = rhombus

        return closest

    def _rhombii(self, distance: float, forward: bool):
        intersection_tuples = []

        lattice_coords = [0]*5
        lattice_coords[self.family.pentangle.pentangle] = self.multiple

        for pentangle in PentAngles.others(self.family.pentangle):
            other = self.family.tiling.strip_family(pentangle)

            initial_intersection = self.intersection_distance_from_point(other.strip(0))
            delta = initial_intersection - distance
            inverse_sin = self.family.pentangle.inverse_sin(pentangle)

            if (forward and inverse_sin > 0) or (not forward and inverse_sin <= 0):
                multiple = int(math.floor(delta / inverse_sin))
                lattice_coords[pentangle.pentangle] = multiple
            else:
                multiple = int(math.ceil(delta / inverse_sin))
                lattice_coords[pentangle.pentangle] = multiple - 1

            intersection = initial_intersection - inverse_sin * multiple

            intersection_tuples.append((pentangle, multiple, intersection))
            intersection_tuples.sort(key=lambda val: val[2], reverse=not forward)

        while True:
            closest_tuple = intersection_tuples.pop(0)

            inverse_sin = self.family.pentangle.inverse_sin(closest_tuple[0])

            if forward and inverse_sin < 0 or not forward and inverse_sin > 0:
                lattice_coords[closest_tuple[0].pentangle] += 1
                rhombus_coords = tuple(lattice_coords)
                new_tuple = (closest_tuple[0], closest_tuple[1]+1, closest_tuple[2] - inverse_sin)
            else:
                rhombus_coords = tuple(lattice_coords)
                lattice_coords[closest_tuple[0].pentangle] -= 1
                new_tuple = (closest_tuple[0], closest_tuple[1]-1, closest_tuple[2] + inverse_sin)

            yield Rhombus(
                self,
                self.family.tiling.strip_family(closest_tuple[0]).strip(closest_tuple[1]),
                rhombus_coords)

            for i in range(2, -1, -1):
                if (forward and new_tuple[2] > intersection_tuples[i][2]) or (
                        not forward and new_tuple[2] < intersection_tuples[i][2]):
                    intersection_tuples.insert(i+1, new_tuple)
                    break

    def origin(self):
        """Returns the origin of the strip, which is an arbitrary point that lies on the strip."""
        return self.family.offset_direction() * (self.family.offset + self.multiple)

    def two_points(self):
        """Returns two arbitrary, non-equal points on the strip as vectors from the origin."""
        point1 = self.origin()
        return point1, point1 + (self.family.direction() * 1000)

    def intersection(self, other: 'Strip'):
        """Returns the coordinates of the intersection between this Strip and the given Strip.

        :param other: The strip to find the intersection with.

        :return: The coordinates in pentagrid space of the intersection between this strip and the given strip, or None
                 if the strips are parallel.
        """
        if self.family.pentangle == other.family.pentangle:
            return None

        return _intersection(*self.two_points(), *other.two_points())

    def intersection_distance_from_point(self, other: 'Strip'):
        """Returns the distance in pentagrid space from the strip's origin to the intersection with the given strip.

        :param other: The strip to find the intersection with. Must not be from the same family/parallel.
        """
        intersection = self.intersection(other)
        if not intersection:
            raise ValueError("The strips don't intersect")

        return (intersection - self.origin()).dot(self.family.direction())

    def __eq__(self, other):
        if not isinstance(other, Strip):
            return False

        return self.family.pentangle == other.family.pentangle and self.multiple == other.multiple

    def __hash__(self):
        return hash((self.family.pentangle, self.multiple))

    def __repr__(self):
        return "Strip(%d, %f)" % (self.family.pentangle.pentangle, self.multiple)


class RhombusVertex(object):
    """Represents a single rhombus vertex."""
    def __init__(self, coordinate, lattice_coordinate):
        """Constructs a new RhombusVertex.

        :param coordinate: The coordinate of the vertex in rhombus space.
        :param lattice_coordinate: The lattice coordinate of the vertex.
        """
        self._coordinate = coordinate
        self._lattice_coordinate = lattice_coordinate

    @property
    def coordinate(self):
        """The coordinate of the vertex in rhombus space."""
        return self._coordinate

    @property
    def lattice_coordinate(self):
        """The lattice coordinate of the vertex."""
        return self._lattice_coordinate

    def __repr__(self):
        return "RhombusVertex(%s, %s)" % (self.coordinate, self.lattice_coordinate)


class Rhombus(object):
    """Represents a single Rhombus."""
    # a list of lattice coordinate offsets, that produce the list of vertices in order around the rhombus
    _VERTEX_OFFSETS = [
        (0, 0),
        (0, -1),
        (-1, -1),
        (-1, 0)]

    def __init__(self, strip1: Strip, strip2: Strip, lattice_coords: Tuple[int, ...]):
        """Constructs a new Rhombus.

        Note that the order of the two strips does not matter.

        :param strip1: One of the pair of intersecting strips that defines this rhombus.
        :param strip2: The other of the pair of intersecting strips that defines this rhombus.
        :param lattice_coords: The lattice coordinates of the "maximal" vertex.
        """
        self._strip1 = strip1
        self._strip2 = strip2
        self._lattice_coords = lattice_coords

    @property
    def strip1(self):
        """One of the pair of intersecting strips that defines this rhombus."""
        return self._strip1

    @property
    def strip2(self):
        """The other of the pair of intersecting strips that defines this rhombus."""
        return self._strip2

    @property
    def lattice_coords(self):
        """The lattice coordinate of the "maximal" vertex."""
        return self._lattice_coords

    def type(self):
        """Returns the RhombusType of this rhombus."""
        diff = abs(self.strip1.family.pentangle.pentangle - self.strip2.family.pentangle.pentangle)

        if diff == 1 or diff == 4:
            return RhombusType.THICK
        return RhombusType.THIN

    def vertices(self):
        """Returns the 4 vertices of this rhombus."""
        vertices = []

        for i in range(0, 4):
            vertex_offsets = Rhombus._VERTEX_OFFSETS[i]

            coords = [*self.lattice_coords]
            coords[self.strip1.family.pentangle.pentangle] += vertex_offsets[0]
            coords[self.strip2.family.pentangle.pentangle] += vertex_offsets[1]

            vertices.append(
                RhombusVertex(self._cartesian_from_lattice(coords), tuple(coords)))

        if _ccw(vertices[0].coordinate, vertices[1].coordinate, vertices[2].coordinate) > 0:
            vertices.reverse()
        return vertices

    @cached_property
    def midpoint(self):
        """Returns the midpoint of this rhombus."""
        vertices = []
        # get 2 opposite vertices
        for i in (0, 2):
            vertex_offsets = Rhombus._VERTEX_OFFSETS[i]

            coords = [*self.lattice_coords]
            coords[self.strip1.family.pentangle.pentangle] += vertex_offsets[0]
            coords[self.strip2.family.pentangle.pentangle] += vertex_offsets[1]

            vertices.append(self._cartesian_from_lattice(coords))
        return Vector(
            (vertices[0].x + vertices[1].x) / 2,
            (vertices[0].y + vertices[1].y) / 2)

    def ordered_strips(self):
        """Returns the 2 intersection Strips that define this rhombus, in order of pentangle."""
        if self.strip1.family.pentangle.pentangle < self.strip2.family.pentangle.pentangle:
            return self.strip1, self.strip2
        return self.strip2, self.strip1

    def contains_point(self, point: Vector):
        """Returns whether this rhombus contains the given point.

        Will return true for any point on the interior, any point on one of the edges, or any of the vertices."""
        vertices = self.vertices()
        for i in range(0, 4):
            val = _ccw(vertices[i].coordinate, vertices[(i+1) % 4].coordinate, point)
            if val > 1e-10:
                return False
        return True

    @staticmethod
    def _cartesian_from_lattice(lattice_coords: List[int]):
        x = 0.0
        y = 0.0

        for pentangle in PentAngles.all():
            x += lattice_coords[pentangle.pentangle] * pentangle.cos()
            y -= lattice_coords[pentangle.pentangle] * pentangle.sin()

        return Vector(x, y)

    def __eq__(self, other):
        if not isinstance(other, Rhombus):
            return False

        return self.ordered_strips() == other.ordered_strips()

    def __hash__(self):
        return hash(self.ordered_strips())

    def __repr__(self):
        return "Rhombus(%s, %s)" % (self.strip1, self.strip2)


class Tiling(object):
    """Represents a P3 penrose tiling, generated by de Bruijn's method

    (1) https://www.mathpages.com/home/kmath621/kmath621.htm
    (2) http://www.ams.org/publicoutreach/feature-column/fcarc-ribbons
    """

    def __init__(self,
                 offsets: Optional[Sequence[float]] = None,
                 rnd: Optional[random.Random] = None,
                 check_offsets: bool = True):
        """Construct a new tiling.

        :param offsets: A list of 5 offsets, as the offset for each strip family.

        In order to generate a proper penrose tiling, the offsets must follow the constraints outlined in (1). They
        must be distinct from each other, the sum of any two must not be an integer, and the sum of all 5 must be an
        integer.

        If not provided, offsets will be generated randomly using rnd.
        :param rnd: A random number generator to use to generate the offsets. Only used if the offsets are not provided.
        If neither offsets nor rnd is not provided, a new Random instance will be created with the default seeding
        behavior.
        :param check_offsets: If true, and the offsets are provided, the offsets will be checked to ensure they meet
        the constraints that will generate a proper penrose tiling.

        If false, offsets not meeting these constraints can be used. This will still generate a valid tiling, and can
        be interesting in its own right, it just won't be a proper P3 penrose tiling.
        """
        if offsets and check_offsets:
            offset_sum = 0
            for i in range(0, 5):
                offset_sum += offsets[i]
                for j in range(i+1, 5):
                    if (offsets[i] + offsets[j]) == int(offsets[i] + offsets[j]):
                        raise ValueError("The sum of offset %d and %d are an integer: %f + %f = %d" % (
                            i, j, offsets[i], offsets[j], int(offsets[i] + offsets[j])))
            if int(offset_sum) != offset_sum:
                raise ValueError("The sum of all 5 offsets must be an integer.")

        if not offsets:
            if not rnd:
                rnd = random.Random()
            offsets = self._generate_offsets(rnd)

        self._families = [StripFamily(self, offsets[pentangle.pentangle], pentangle) for pentangle in PentAngles.all()]

    @staticmethod
    def _generate_offsets(rnd: random.Random):
        offsets = []

        # it's theoretically possible that 2 of the randomly generated values below will have a sum that is an
        # integer, but in practice it shouldn't be an issue
        offset_sum = 0
        for i in range(0, 4):
            offsets.append(rnd.random())
            offset_sum += offsets[-1]

        # generate the last offset such that the sum of all 5 offsets is an integer, per the constraints
        offsets.append(math.ceil(offset_sum) - offset_sum)

        return offsets

    def strip_family(self, pentangle: PentAngle):
        """Returns the StripFamily for the given pentangle.

        :param pentangle: The pentangle of the StripFamily.
        """
        return self._families[pentangle.pentangle]

    def rhombus_at_point(self, point: Vector):
        """Returns the Rhombus that contains the given point.

        :param point: The point to find the containing Rhombus of.
        """
        strips: List[Strip] = []
        # find all strips that are near enough to the point that could feasibly contain the point
        for pentangle in PentAngles.all():
            family = self.strip_family(pentangle)
            strips.extend(family.strips_near_point(point))

        # iterate over all the intersections between the potential strips
        for strip, other_strip in itertools.combinations(strips, 2):
            if strip.family == other_strip.family:
                continue

            rhombus = strip.rhombus_at_intersection(other_strip)
            if rhombus.contains_point(point):
                return rhombus
        raise Exception("Could not find the rhombus containing the given point")

    def rhombus_edges(self, grid_cell: GridCell):
        """Returns all distinct rhombus edges of rhombii in the given grid cell.

        All edges will be returned for every rhombus whose midpoint is contained in the grid cell. Some of the returned
        edges will fall outside the grid cell.

        A rhombus is contained by exactly one grid cell. If the midpoint lies on a grid cell edge or corner, the grid
        cell with the higher x and/or y indices "wins".

        :param grid_cell: The grid cell to find the rhombus edges within.
        :returns: An iterable of (RhombusVertex, RhombusVertex) tuples.
        """
        processed_edges = set()

        for rhombus in self.rhombii(grid_cell):
            vertices = rhombus.vertices()
            for i in range(4):
                edge_key = (
                    vertices[i].lattice_coordinate,
                    vertices[(i+1) % 4].lattice_coordinate)

                if edge_key not in processed_edges and tuple(reversed(edge_key)) not in processed_edges:
                    processed_edges.add(edge_key)
                    edge = vertices[i], vertices[(i+1) % 4]
                    yield edge

    def rhombii(self, grid_cell: GridCell):
        """Returns all rhombii in the given grid cell.

        All rhombii whose midpoint lies in the given grid cell will be returned. Parts of some of these rhombii will
        extend outside the grid cell.

        A rhombus is contained by exactly one grid cell. If the midpoint lies on a grid cell edge or corner, the grid
        cell with the higher x and/or y indices "wins".

        :param grid_cell: The grid cell to get the rhombii within.
        """

        def rhombus_in_cell(r: Rhombus):
            midpoint = r.midpoint
            return (grid_cell.origin.x <= midpoint.x < grid_cell.extent.x and
                    grid_cell.origin.y <= midpoint.y < grid_cell.extent.y)

        grid_corners = [*grid_cell.corners(1.6)]
        for pentangle in PentAngles.others(PentAngles.pentangle(4)):
            family = self.strip_family(pentangle)

            min_multiple = math.inf
            max_multiple = -math.inf
            for corner in grid_corners:
                for strip in family.strips_near_point(corner):
                    if strip.multiple < min_multiple:
                        min_multiple = strip.multiple
                    if strip.multiple > max_multiple:
                        max_multiple = strip.multiple

            for multiple in range(min_multiple, max_multiple + 1):
                strip = family.strip(multiple)

                intersection_distances = []
                two_points = strip.two_points()
                for corner1, corner2 in itertools.pairwise(grid_corners + [grid_corners[0]]):
                    if _ccw(*two_points, corner1/2.5) != _ccw(*two_points, corner2/2.5):
                        intersection = _intersection(*two_points, corner1/2.5, corner2/2.5)
                        intersection_distances.append(
                            (intersection - strip.origin()).dot(strip.family.direction()))

                if len(intersection_distances) != 2:
                    continue

                if intersection_distances[0] > intersection_distances[1]:
                    intersection_distances.reverse()

                start_rhombus = strip.rhombus(intersection_distances[0])
                stop_rhombus = strip.rhombus(intersection_distances[1])
                if start_rhombus == stop_rhombus:
                    if rhombus_in_cell(start_rhombus):
                        yield start_rhombus
                    continue

                for rhombus in strip.rhombii(intersection_distances[0], True):
                    # since we're processing the strips in order by family, we can ignore any rhombus from
                    # an intersection with a previous family
                    if rhombus.strip2.family.pentangle.pentangle > strip.family.pentangle.pentangle:
                        if rhombus_in_cell(rhombus):
                            yield rhombus

                    if rhombus == stop_rhombus:
                        break
