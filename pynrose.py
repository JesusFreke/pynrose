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
from functools import total_ordering, cached_property
from typing import Iterator, Iterable, List, Optional


class Vector(object):
    def __init__(self, x, y):
        self._x = x
        self._y = y

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @cached_property
    def length(self):
        return math.sqrt(self._x ** 2 + self._y ** 2)

    def dot(self, other: 'Vector'):
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
    def __init__(self, origin: Vector, grid_size: Vector):
        self.origin = origin
        self.grid_size = grid_size

    def cell(self, x: int, y: int):
        return GridCell(self, x, y)

    def __eq__(self, other):
        if not isinstance(other, Grid):
            return False
        return self.origin == other.origin and self.grid_size == other.grid_size

    def __hash__(self):
        return hash((self.origin, self.grid_size))


@total_ordering
class GridCell(object):
    def __init__(self, grid: Grid, x_multiple: int, y_multiple: int):
        self.grid = grid
        self.x_multiple = x_multiple
        self.y_multiple = y_multiple

        self.origin = grid.origin + Vector(grid.grid_size.x * x_multiple, grid.grid_size.y * y_multiple)
        self.extent = self.origin + grid.grid_size

    def midpoint(self):
        return (self.origin + self.extent) / 2

    def __eq__(self, other):
        if not isinstance(other, GridCell):
            return False

        return (self.grid, self.x, self.y) == (other.grid, other.x, other.y)

    def __lt__(self, other):
        if not isinstance(other, GridCell):
            raise TypeError("Cannot compare GridCell with %s" % type(other))

        if self.origin.x < other.origin.x:
            return True
        elif self.origin.x > other.origin.x:
            return False

        if self.origin.y < other.origin.y:
            return True

        return False

    def __hash__(self):
        return hash((self.grid, self.x, self.y))


class RhombusType(Enum):
    THIN = 0
    THICK = 1


class PentAngle(object):
    def __init__(self, pentangle: int):
        self.pentangle: int = int(pentangle % 5)

        self.radians = self.pentangle * 2 * math.pi / 5

        self._sin = math.sin(self.radians)
        self._cos = math.cos(self.radians)

    def unit(self) -> 'Vector':
        """Returns a unit vector in the direction of this PentAngle"""
        return Vector(self.sin(), self.cos())

    def sin(self, other: 'PentAngle' = None) -> float:
        """Returns the sin of the angle between this pentangle and another pentangle"""
        if other is None:
            return self._sin

        return math.sin(other.radians - self.radians)

    def cos(self, other: 'PentAngle' = None) -> float:
        """Returns the cos of the angle between this pentangle and another pentangle"""
        if other is None:
            return self._cos

        return math.cos(other.radians - self.radians)

    def __eq__(self, other):
        if not isinstance(other, PentAngle):
            return False
        return other.pentangle == self.pentangle

    def __hash__(self):
        return self.pentangle

    @staticmethod
    def all() -> Iterable['PentAngle']:
        return [
            PentAngle(0),
            PentAngle(1),
            PentAngle(2),
            PentAngle(3),
            PentAngle(4)]

    @staticmethod
    def others(pent_angle: 'PentAngle') -> Iterator['PentAngle']:
        """Yields the other 4 pentangles, excluding the given one."""
        for other in PentAngle.all():
            if other == pent_angle:
                continue
            yield other


class StripFamily(object):
    """Represents 1 of the 5 families of strips of rhombii in de Bruijn's method"""

    def __init__(self, tiling: 'Tiling', offset: float, pentangle: PentAngle):
        self.tiling = tiling
        self.offset = offset
        self.pentangle = pentangle

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
        return Strip(self, multiple)

    def origin(self):
        return self.strip(0).origin()

    def strips_near_point(self, point: Vector):
        """Returns 1 or 2 strips from this family that are near the given point.

        The returned strips are ones that could potentially contain the rhombus at the given point (but may not).
        """
        pentagrid_point = point / 2.5
        multiple = (pentagrid_point - self.origin()).dot(self.offset_direction())

        # the rhombii associated with a strip may wander away from the strip up to a max of just under 2 units in
        # rhombus space, which is .8 units in pentagrid space
        if abs(math.ceil(multiple) - multiple) <= .8:
            yield self.strip(math.ceil(multiple))
        if abs(math.floor(multiple) - multiple) <= .8:
            yield self.strip(math.floor(multiple))

    def strip_range_containing_points(self, points: List[Vector]):
        pass


def _det(x1, y1, x2, y2):
    """Returns the determinate of (x1, y1) and (x2, y2)

    | x1 y1 |
    | x2 y2 |
    """
    return x1 * y2 - y1 * x2


def _ccw(p1, p2, p3):
    return (p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y)


class Strip(object):
    """Represents a single string of rhombii in de Bruijn's method"""
    def __init__(self, family: StripFamily, multiple: int):
        self.family = family
        self.multiple = multiple

    def rhombii(self, start: 'Strip', forward: bool):
        """Returns a generator over the rhombii for this strip, starting at the intersection with 'start'.

        :param start: A non-parallel strip. The rhombus representing the intersection of this strip with the given
        start strip will be used as the starting point of the iteration
        :param forward: The direction to iterate
        :return: An infinite generator starting at the intersection of this strip and the given starting strip, and in
        the given direction
        """
        return self._rhombii(start, self.intersection_distance_from_point(start), forward)

    def rhombus(self, target_distance: float):
        """Gets the rhombus near the given distance from this strip's origin.

        This finds the intersection closest to the point that is the given distance from this strip's origin, and
        returns the associated rhombus
        """
        closest_value = math.inf
        closest = None

        # start iterating just a bit before the target distance
        rhombii = self._rhombii(None, target_distance, True)

        # TODO: do some logging and experiments to see how far this actually needs to go
        # 10 should be plenty to ensure we find the closest intersections on either side of the target point
        for rhombus in [next(rhombii) for i in range(10)]:
            other_strip = rhombus.strip2 if rhombus.strip1 == self else rhombus.strip1
            distance_from_target = abs(self.intersection_distance_from_point(other_strip) - target_distance)
            if distance_from_target < closest_value:
                closest_value = distance_from_target
                closest = rhombus
        return closest

    def _rhombii(self, start: Optional['Strip'], distance: float, forward: bool):
        intersections = [0.0] * 5
        intersection_multiples = [0] * 5

        for pentangle in PentAngle.others(self.family.pentangle):
            if start is not None and pentangle == start.family.pentangle:
                intersections[pentangle.pentangle] = distance
                intersection_multiples[pentangle.pentangle] = start.multiple
                continue

            other = self.family.tiling.strip_family(pentangle)

            initial_intersection = self.intersection_distance_from_point(other.strip(0))
            delta = initial_intersection - distance
            sin = pentangle.sin(self.family.pentangle)
            interval = 1/sin

            # TODO: do we need to do the split ceil/floor thing? can we always use one or the other?
            # also, why negative?
            if forward:
                if interval < 0:
                    intersection_multiples[pentangle.pentangle] = -int(math.ceil(delta / interval))
                else:
                    intersection_multiples[pentangle.pentangle] = -int(math.floor(delta / interval))
            else:
                if interval < 0:
                    intersection_multiples[pentangle.pentangle] = -int(math.floor(delta / interval))
                else:
                    intersection_multiples[pentangle.pentangle] = -int(math.ceil(delta / interval))

            intersections[pentangle.pentangle] = (
                self.intersection_distance_from_point(other.strip(intersection_multiples[pentangle.pentangle])))

        while True:
            closest = None
            closest_value = math.inf if forward else -math.inf

            for pentangle in PentAngle.others(self.family.pentangle):
                if ((forward and intersections[pentangle.pentangle] < closest_value) or
                        (not forward and intersections[pentangle.pentangle] > closest_value)):
                    closest_value = intersections[pentangle.pentangle]
                    closest = pentangle

            lattice_coords = [0] * 5

            for pentangle in PentAngle.all():
                if pentangle == self.family.pentangle:
                    lattice_coords[pentangle.pentangle] = self.multiple
                elif pentangle == closest:
                    lattice_coords[pentangle.pentangle] = intersection_multiples[pentangle.pentangle]
                else:
                    lattice_coords[pentangle.pentangle] = intersection_multiples[pentangle.pentangle]

                    sin = self.family.pentangle.sin(pentangle)
                    if (forward and sin < 0) or (not forward and sin > 0):
                        lattice_coords[pentangle.pentangle] -= 1

            closest_multiple = intersection_multiples[closest.pentangle]
            sin = closest.sin(self.family.pentangle)

            if (forward and sin < 0) or (not forward and sin > 0):
                intersection_multiples[closest.pentangle] -= 1
            else:
                intersection_multiples[closest.pentangle] += 1

            intersections[closest.pentangle] = self.intersection_distance_from_point(
                self.family.tiling.strip_family(closest).strip(intersection_multiples[closest.pentangle]))

            yield Rhombus(
                self,
                self.family.tiling.strip_family(closest).strip(closest_multiple),
                lattice_coords)

    def origin(self):
        """Gets the origin of the strip, which is an arbitrary point that lies on the strip"""
        return self.family.offset_direction() * (self.family.offset + self.multiple)

    def two_points(self):
        """Gets two arbitrary, non-equal points on the strip as vectors from the origin"""
        point1 = self.origin()
        return point1, point1 + (self.family.direction() * 1000)

    def intersection(self, other: 'Strip'):
        """Gets the coordinates of the intersection between this strip and the given strip, or None if parallel"""

        if self.family.pentangle == other.family.pentangle:
            return None

        # see, e.g. https://mathworld.wolfram.com/Line-LineIntersection.html

        (p0, p1) = self.two_points()
        (p2, p3) = other.two_points()

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

        return Vector(x_num / denom, y_num / denom)

    def intersection_distance_from_point(self, other: 'Strip'):
        """Gets the distance from the strip's origin to the intersection with the given strip"""
        intersection = self.intersection(other)
        if not intersection:
            raise ValueError("The strips don't intersect")

        # TODO: probably need to check if some small distance from 0, because floating point
        if intersection.length == 0:
            return 0
        
        return (intersection - self.origin()).dot(self.family.direction())

    def __eq__(self, other):
        if not isinstance(other, Strip):
            return False

        return self.family.pentangle == other.family.pentangle and self.multiple == other.multiple

    def __hash__(self):
        return hash((self.family.pentangle, self.multiple))

    def __repr__(self):
        return "Strip(%d, %f)" % (self.family.pentangle.pentangle, self.multiple)


class Rhombus(object):
    # a list of lattice coordinate offsets, that produce the list of vertices in the correct order around the rhombus
    _VERTEX_OFFSETS = [
        (0, 0),
        (0, -1),
        (-1, -1),
        (-1, 0)]

    def __init__(self, strip1: Strip, strip2: Strip, lattice_coords: List[int]):
        self.strip1 = strip1
        self.strip2 = strip2
        self.lattice_coords = lattice_coords

    def type(self):
        diff = abs(self.strip1.family.pentangle.pentangle - self.strip2.family.pentangle.pentangle)

        if diff == 1 or diff == 4:
            return RhombusType.THICK
        return RhombusType.THIN

    def vertices(self):
        vertices = []

        for i in range(0, 4):
            vertex_offsets = Rhombus._VERTEX_OFFSETS[i]

            coords = self.lattice_coords.copy()
            coords[self.strip1.family.pentangle.pentangle] += vertex_offsets[0]
            coords[self.strip2.family.pentangle.pentangle] += vertex_offsets[1]

            vertices.append(self._cartesian_from_lattice(coords))

        if _ccw(vertices[0], vertices[1], vertices[2]) > 0:
            vertices.reverse()
        return vertices

    @cached_property
    def midpoint(self):
        vertices = []
        # get 2 opposite vertices
        for i in (0, 2):
            vertex_offsets = Rhombus._VERTEX_OFFSETS[i]

            coords = self.lattice_coords.copy()
            coords[self.strip1.family.pentangle.pentangle] += vertex_offsets[0]
            coords[self.strip2.family.pentangle.pentangle] += vertex_offsets[1]

            vertices.append(self._cartesian_from_lattice(coords))
        return Vector(
            (vertices[0].x + vertices[1].x) / 2,
            (vertices[0].y + vertices[1].y) / 2)

    def ordered_strips(self):
        """returns the 2 intersection Strips that define this Rhombus, in order of pentangle"""
        if self.strip1.family.pentangle.pentangle < self.strip2.family.pentangle.pentangle:
            return self.strip1, self.strip2
        return self.strip2, self.strip1

    def contains_point(self, point: Vector):
        vertices = self.vertices()
        for i in range(0, 4):
            val = _ccw(vertices[i], vertices[(i+1) % 4], point)
            if val > 1e-10:
                return False
        return True

    @staticmethod
    def _cartesian_from_lattice(lattice_coords: List[int]):
        x = 0.0
        y = 0.0

        for pentangle in PentAngle.all():
            x += lattice_coords[pentangle.pentangle] * pentangle.cos()
            y -= lattice_coords[pentangle.pentangle] * pentangle.sin()

        return Vector(x, y)

    def __eq__(self, other):
        if not isinstance(other, Rhombus):
            return False

        return self.ordered_strips() == other.ordered_strips()

    def __hash__(self):
        return hash(self.ordered_strips())


class Tiling(object):
    """Represents a P3 penrose tiling, generated by de Bruijn's method

    (1) https://www.mathpages.com/home/kmath621/kmath621.htm
    (2) http://www.ams.org/publicoutreach/feature-column/fcarc-ribbons
    """

    def __init__(self, offsets: Optional[List[float]]=None, rnd: Optional[random.Random]=None):
        """Create a new tiling.

        :param offsets: A list of 5 offsets, as the offset for each strip family. These offsets must follow the
            constraints outlined in (1). They must be distinct from each other, the sum of any two must not be an
            integer, and the sum of all 5 must be an integer.

            If not provided, offsets will be generated randomly using rnd.
        :param rnd: A random number generator to use to generate the offsets. Only used if offsets is None. If not
        provided, a new Random instance will be created with the default seeding behavior.
        """
        # TODO: check the offset constraints
        # TODO: normalize the offsets to [0, 1)
        if not offsets:
            if not rnd:
                rnd = random.Random()
            offsets = self._generate_offsets(rnd)

        self._families = [StripFamily(self, offsets[pentangle.pentangle], pentangle) for pentangle in PentAngle.all()]

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
        return self._families[pentangle.pentangle]

    def rhombus_at_point(self, point: Vector):
        """Returns the rhombus that contains the given point"""
        pentagrid_point = point/2.5

        strips: List[Strip] = []
        # find all strips that are near enough to the point that could feasibly contain the point
        for pentangle in PentAngle.all():
            family = self.strip_family(pentangle)
            strips.extend(family.strips_near_point(point))

        # iterate over all the intersections between the potential strips
        for strip, other_strip in itertools.combinations(strips, 2):
            if strip.family == other_strip.family:
                continue

            rhombus = next(strip.rhombii(other_strip, True))
            if rhombus.contains_point(point):
                return rhombus
        raise Exception("Could not find the rhombus containing the given point")

    def rhombii(self, grid_cell: GridCell):
        processed_rhombii = set()
        pending_rhombii = []
        processed_strips = set()

        cell_midpoint = grid_cell.midpoint()
        # every unit in pentagrid space is ~2.5 units in the penrose space
        approximate_multiple = int(cell_midpoint.x / 2.5)
        initial_strip = self.strip_family(PentAngle(0)).strip(approximate_multiple)
        approximate_distance = ((cell_midpoint/2.5) - initial_strip.origin()).dot(initial_strip.family.direction())
        initial_rhombus = initial_strip.rhombus(approximate_distance)

        def rhombus_in_cell(rhombus: Rhombus):
            midpoint = rhombus.midpoint
            if midpoint.x < grid_cell.origin.x or midpoint.x >= grid_cell.extent.x:
                return False
            if midpoint.y < grid_cell.origin.y or midpoint.y >= grid_cell.extent.y:
                return False
            return True

        if not rhombus_in_cell(initial_rhombus):
            raise Exception("The initial rhombus is outside the bounding box. Maybe the bounding boxes are too small?")

        def process_rhombus(rhombus: Rhombus):
            if rhombus in processed_rhombii:
                # TODO: in the original code, we added the rhombus to pending_rhombii regardless. is this needed?
                return
            processed_rhombii.add(rhombus)
            yield rhombus
            pending_rhombii.append(rhombus)

        def continue_processing_strip(rhombus: Rhombus):
            # +/- 2.5, in order to catch the case of a strip parallel with an edge that goes in and out of the grid cell.
            midpoint = rhombus.midpoint
            if (midpoint.x < (grid_cell.origin.x - 2.5)) or (midpoint.x > (grid_cell.extent.x + 2.5)):
                return False
            if (midpoint.y < (grid_cell.origin.y - 2.5)) or (midpoint.y > (grid_cell.extent.y + 2.5)):
                return False
            return True

        def process_strip(strip1: Strip, strip2: Strip):
            if strip1 in processed_strips:
                return
            processed_strips.add(strip1)

            first = True
            for rhombus in strip1.rhombii(strip2, True):
                if first:
                    first = False
                    continue
                if rhombus_in_cell(rhombus):
                    yield from process_rhombus(rhombus)
                else:
                    processed_rhombii.add(rhombus)

                if not continue_processing_strip(rhombus):
                    break

            first = True
            for rhombus in strip1.rhombii(strip2, False):
                if first:
                    first = False
                    continue
                if rhombus_in_cell(rhombus):
                    yield from process_rhombus(rhombus)
                else:
                    processed_rhombii.add(rhombus)

                if not continue_processing_strip(rhombus):
                    break

        yield from process_rhombus(initial_rhombus)
        yield from process_strip(initial_rhombus.strip1, initial_rhombus.strip2)
        yield from process_strip(initial_rhombus.strip2, initial_rhombus.strip1)

        while pending_rhombii:
            pending_rhombus = pending_rhombii.pop(0)
            yield from process_strip(pending_rhombus.strip2, pending_rhombus.strip1)


def print_rhombus_svg(rhombus: Rhombus):
    string = '<path'

    if rhombus.type() == RhombusType.THIN:
        string += ' class="thinRhombus"'
    else:
        string += ' class="thickRhombus"'

    string += ' id="Rhombus (%d, %d) (%d, %d)"' % (
        rhombus.strip1.family.pentangle.pentangle,
        rhombus.strip1.multiple,
        rhombus.strip2.family.pentangle.pentangle,
        rhombus.strip2.multiple)

    string += ' d="M'

    for vertex in rhombus.vertices():
        string += ' %f,%f' % (vertex.x, vertex.y)
    string += ' z"/>'
    print(string)


def main():
    tiling = Tiling(rnd=random.Random(12345))

    grid = Grid(Vector(-100, -100), Vector(100, 100))

    max_protrusion = math.sin(math.radians(72))

    view_size = Vector(
        200 + max_protrusion * 2,
        200 + max_protrusion * 2)

    print('<svg width="%fmm" height="%fmm" viewBox="%f %f %f %f">' % (
        view_size.x,
        view_size.y,
        -view_size.x/2,
        -view_size.y/2,
        view_size.x/2 + max_protrusion,
        view_size.y/2 + max_protrusion))
    print('<style><![CDATA[')
    print('rect.boundingBox {')
    print('    stroke: blue;')
    print('    stroke-width: .05;')
    print('    fill-opacity: 0;')
    print('    stroke-opacity: .5;')
    print('}')
    print('path.thinRhombus {')
    print('    fill: #333333;')
    print('    stroke: #000000;')
    print('    stroke-width: .01;')
    print('}')
    print('path.thickRhombus {')
    print('    fill: #aaaaaa;')
    print('    stroke: #000000;')
    print('    stroke-width: .01;')
    print('}')
    print(']]></style>')

    for rhombus in tiling.rhombii(grid.cell(0, 0)):
        print_rhombus_svg(rhombus)

    print('<rect x="%f" y="%f" width="%f" height="%f" class="boundingBox"/>' % (
        -100, -100, 100, 100))

    print('</svg>')


if __name__ == "__main__":
    main()






        
