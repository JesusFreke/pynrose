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
import multiprocessing
import random
import unittest
from multiprocessing import Process
from typing import List

from pynrose import Tiling, Vector, Grid, PentAngles, Rhombus, GridCell


class TestTiling(unittest.TestCase):

    def test_rhombus_at_intersection(self):
        tiling = Tiling(rnd=random.Random(12345))

        for rhombus in tiling.rhombii(Grid(Vector(0, 0), Vector(100, 100)).cell(0, 0)):
            intersection_rhombus = rhombus.strip1.rhombus_at_intersection(rhombus.strip2)
            assert intersection_rhombus == rhombus
            assert intersection_rhombus.lattice_coords == rhombus.lattice_coords

    def test_sins(self):
        fifth_angle = math.pi * 2 / 5

        for p1, p2 in itertools.permutations(PentAngles.all(), 2):
            assert math.isclose(p1.sin(p2), math.sin((p2.pentangle - p1.pentangle) * fifth_angle))

    def test_rhombus_contains_point(self):
        tiling = Tiling(rnd=random.Random(12345))
        rhombus = next(tiling.rhombii(Grid(Vector(0, 0), Vector(20, 20)).cell(0, 0)))

        vertices = rhombus.vertices()

        assert rhombus.contains_point(rhombus.midpoint)
        for vertex in rhombus.vertices():
            assert rhombus.contains_point(vertex)

        assert not rhombus.contains_point(Vector(rhombus.midpoint.x + 3, rhombus.midpoint.y))
        assert not rhombus.contains_point(Vector(rhombus.midpoint.x - 3, rhombus.midpoint.y))
        assert not rhombus.contains_point(Vector(rhombus.midpoint.x, rhombus.midpoint.y + 3))
        assert not rhombus.contains_point(Vector(rhombus.midpoint.x, rhombus.midpoint.y - 3))

        # a point just barely past a vertex
        assert not rhombus.contains_point(rhombus.midpoint + (vertices[0] - rhombus.midpoint) * 1.01)

        midpoint = vertices[0] + (vertices[1] - vertices[0]) / 2
        assert rhombus.contains_point(midpoint)

        # a point just barely past the midpoint of an edge
        assert not rhombus.contains_point(rhombus.midpoint + (midpoint - rhombus.midpoint) * 1.01)

        # a point just barely past the vertex
        assert not rhombus.contains_point(vertices[0] + (vertices[0] - rhombus.midpoint) * 1.01)

    def test_rhombus_at_point(self):
        tiling = Tiling(rnd=random.Random(12345))

        rnd = random.Random(12345)

        for i in range(100):
            target = Vector(rnd.uniform(-1000, 1000), rnd.uniform(-1000, 1000))

            rhombus = tiling.rhombus_at_point(target)
            if not rhombus.contains_point(target):
                self.fail("%s is not in %s" % (target, rhombus.vertices()))

            assert tiling.rhombus_at_point(rhombus.midpoint) == rhombus

            vertices = rhombus.vertices()
            edge_midpoint = vertices[0] + (vertices[0] - vertices[1]) / 2

            edge_rhombus = tiling.rhombus_at_point(edge_midpoint)
            assert edge_rhombus.contains_point(edge_midpoint)

            for vertex in rhombus.vertices():
                rhombus = tiling.rhombus_at_point(vertex)
                assert rhombus.contains_point(vertex)

    def test_backward_forward_strip_iteration(self):
        tiling = Tiling(rnd=random.Random(12345))

        rnd = random.Random(12345)

        for pentangle in PentAngles.all():
            for _ in range(100):
                strip = tiling.strip_family(pentangle).strip(rnd.randint(-100, 100))

                rhombii: List[Rhombus] = [*itertools.islice(
                    strip.rhombii(
                        rnd.uniform(-100, 100),
                        True),
                    20)]

                backwards_rhombii = [*itertools.islice(
                    strip.rhombii(
                        strip.intersection_distance_from_point(rhombii[-1].strip2) + 1e-9,
                        False),
                    20)]

                backwards_rhombii.reverse()

                assert rhombii == backwards_rhombii

                # check that each rhombus shares 2 vertices with the next rhombus in the strip
                for rhombus1, rhombus2 in itertools.pairwise(rhombii):
                    shared_vertex_count = 0
                    for rhombus1_vertex in rhombus1.vertices():
                        for rhombus2_vertex in rhombus2.vertices():
                            if (math.isclose(rhombus1_vertex.x, rhombus2_vertex.x) and
                                    math.isclose(rhombus1_vertex.y, rhombus2_vertex.y)):
                                shared_vertex_count += 1
                    assert shared_vertex_count == 2

    def test_tiling(self):
        tiling = Tiling(rnd=random.Random(12345))

        tiling_rhombii = [*tiling.rhombii(Grid(Vector(0, 0), Vector(5, 5)).cell(0, 0))]

        expected_rhombii = {((2, -3), (3, 0), (0, -2, -3, 0, 1)), ((1, -1), (2, -3), (1, -1, -3, -1, 1)),
                            ((0, 1), (1, 0), (1, 0, -3, -2, 0)), ((1, -2), (2, -2), (-1, -2, -2, 0, 1)),
                            ((0, 1), (4, 1), (1, -1, -3, -1, 1)), ((0, 0), (4, 0), (0, -1, -2, -1, 0)),
                            ((0, 1), (1, -1), (1, -1, -3, -1, 1)), ((1, -1), (2, -2), (0, -1, -2, -1, 0)),
                            ((2, -1), (3, 0), (-1, -1, -1, 0, -1)), ((2, -2), (3, 0), (0, -2, -2, 0, 0)),
                            ((0, 0), (3, 0), (0, -1, -2, 0, 0)), ((0, 0), (1, -1), (0, -1, -2, 0, 0)),
                            ((2, -2), (3, -1), (0, -1, -2, -1, 0)), ((3, 0), (4, 0), (-1, -1, -2, 0, 0)),
                            ((2, -2), (4, 1), (0, -2, -2, 0, 1)), ((1, 0), (2, -2), (0, 0, -2, -2, 0)),
                            ((1, -1), (4, 1), (0, -1, -3, -1, 1)), ((1, -1), (3, 0), (0, -1, -2, 0, 0)),
                            ((1, 0), (4, 0), (0, 0, -2, -1, 0)), ((1, 0), (3, -1), (0, 0, -2, -1, 0)),
                            ((0, 1), (2, -2), (1, 0, -2, -2, 0)), ((0, 0), (4, 1), (0, -2, -2, 0, 1)),
                            ((0, 0), (2, -2), (0, -2, -2, 0, 1)), ((0, 1), (3, -1), (1, -1, -3, -1, 0)),
                            ((0, 1), (2, -3), (1, -2, -3, -1, 1)), ((3, -1), (4, 1), (1, -1, -3, -1, 1)),
                            ((0, 0), (1, 0), (0, 0, -2, -1, -1)), ((3, 0), (4, 1), (0, -2, -3, 0, 1))}

        actual_rhombii = set()
        for rhombus in tiling_rhombii:
            actual_rhombii.add(((
                    rhombus.strip1.family.pentangle.pentangle,
                    rhombus.strip1.multiple
                ),
                (
                    rhombus.strip2.family.pentangle.pentangle,
                    rhombus.strip2.multiple
                ),
                rhombus.lattice_coords))
            
        assert actual_rhombii == expected_rhombii

        # move the origin just a bit, this should slightly change the rhombii included at the edges
        tiling_rhombii = [*tiling.rhombii(Grid(Vector(.01, .01), Vector(5, 5)).cell(0, 0))]

        expected_rhombii = {((2, -3), (3, 0), (0, -2, -3, 0, 1)), ((1, -1), (2, -3), (1, -1, -3, -1, 1)),
                            ((0, 1), (1, 0), (1, 0, -3, -2, 0)), ((1, -2), (2, -2), (-1, -2, -2, 0, 1)),
                            ((0, 1), (4, 1), (1, -1, -3, -1, 1)), ((0, 0), (4, 0), (0, -1, -2, -1, 0)),
                            ((0, 1), (1, -1), (1, -1, -3, -1, 1)), ((1, -1), (2, -2), (0, -1, -2, -1, 0)),
                            ((2, -2), (3, 0), (0, -2, -2, 0, 0)), ((0, 0), (3, 0), (0, -1, -2, 0, 0)),
                            ((0, 0), (1, -1), (0, -1, -2, 0, 0)), ((2, -2), (3, -1), (0, -1, -2, -1, 0)),
                            ((3, 0), (4, 0), (-1, -1, -2, 0, 0)), ((2, -2), (4, 1), (0, -2, -2, 0, 1)),
                            ((1, 0), (2, -2), (0, 0, -2, -2, 0)), ((1, -1), (4, 1), (0, -1, -3, -1, 1)),
                            ((1, -1), (3, 0), (0, -1, -2, 0, 0)), ((1, 0), (4, 0), (0, 0, -2, -1, 0)),
                            ((1, 0), (3, -1), (0, 0, -2, -1, 0)), ((0, 1), (2, -2), (1, 0, -2, -2, 0)),
                            ((0, 0), (4, 1), (0, -2, -2, 0, 1)), ((0, 0), (2, -2), (0, -2, -2, 0, 1)),
                            ((0, 1), (3, -1), (1, -1, -3, -1, 0)), ((0, 1), (2, -3), (1, -2, -3, -1, 1)),
                            ((3, -1), (4, 1), (1, -1, -3, -1, 1)), ((0, 0), (1, 0), (0, 0, -2, -1, -1)),
                            ((3, 0), (4, 1), (0, -2, -3, 0, 1))}

        actual_rhombii = set()
        for rhombus in tiling_rhombii:
            actual_rhombii.add(((
                                    rhombus.strip1.family.pentangle.pentangle,
                                    rhombus.strip1.multiple
                                ),
                                (
                                    rhombus.strip2.family.pentangle.pentangle,
                                    rhombus.strip2.multiple
                                ),
                                rhombus.lattice_coords))

        assert actual_rhombii == expected_rhombii

    def test_strip_rhombus(self):
        tiling = Tiling(rnd=random.Random(12345))

        rnd = random.Random(12345)

        for pentangle in PentAngles.all():
            for _ in range(100):
                strip = tiling.strip_family(pentangle).strip(rnd.randint(-100, 100))

                rhombii: List[Rhombus] = [*itertools.islice(
                    strip.rhombii(
                        rnd.uniform(-100, 100),
                        True),
                    20)]

                for rhombus1, rhombus2 in itertools.pairwise(rhombii):

                    distance1 = strip.intersection_distance_from_point(rhombus1.strip2)
                    distance2 = strip.intersection_distance_from_point(rhombus2.strip2)

                    assert strip.rhombus(distance1) == rhombus1
                    assert strip.rhombus(distance2) == rhombus2
                    assert strip.rhombus((distance1 + distance2)/2) in (rhombus1, rhombus2)
                    assert strip.rhombus(distance1 + (distance2 - distance1) * .1) == rhombus1
                    assert strip.rhombus(distance1 + (distance2 - distance1) * .9) == rhombus2

    def test_strip_grid_cell_corner_edge_case(self):
        """Tests an edge case where a strip intersects very closely with a grid cell corner, causing an infinite loop"""

        queue = multiprocessing.Queue()

        def run_test():
            tiling = Tiling(rnd=random.Random(1234))
            grid = Grid(Vector(0, 40), Vector(20, 20))
            count = 0
            for _ in tiling.rhombii(grid.cell(0, 0)):
                count += 1
            queue.put(count)

        p = Process(target=run_test)
        p.start()
        p.join(timeout=5)
        if p.exitcode is None:
            p.kill()
            self.fail("rhombus interation encountered an infinite loop")
        rhombus_count = queue.get_nowait()
        assert rhombus_count == 482

    def test_rhombii_at_edges(self):
        """Test that there are no missing rhombii at the edges of a grid cell."""

        tiling = Tiling(rnd=random.Random(12345))
        rnd = random.Random(12345)

        for _ in range(100):
            inner_grid = Grid(
                Vector(rnd.uniform(-100, 100), rnd.uniform(-100, 100)),
                Vector(5, 5))

            # make a grid cell that is bigger than inner_grid by 2.5 on all sides
            outer_grid = Grid(
                inner_grid.origin - Vector(2.5, 2.5),
                inner_grid.grid_size + Vector(5, 5))

            inner_grid_cell = inner_grid.cell(0, 0)
            inner_rhombii = set(tiling.rhombii(inner_grid_cell))

            def rhombus_in_cell(r: Rhombus, cell: GridCell):
                midpoint = r.midpoint
                return (cell.origin.x <= midpoint.x < cell.extent.x and
                        cell.origin.y <= midpoint.y < cell.extent.y)

            for outer_rhombus in tiling.rhombii(outer_grid.cell(0, 0)):
                if rhombus_in_cell(outer_rhombus, inner_grid_cell):
                    assert outer_rhombus in inner_rhombii
