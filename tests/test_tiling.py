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

import random
import unittest
from pynrose import Tiling, Vector, Grid


class TestTiling(unittest.TestCase):

    def test_rhombus_contains_point(self):
        tiling = Tiling(rnd=random.Random(12345))
        rhombus = next(tiling.rhombii(Grid(Vector(0, 0), Vector(20, 20)).cell(0, 0)))

        vertices = rhombus.vertices()

        assert rhombus.contains_point(rhombus.midpoint)
        for vertex in rhombus.vertices():
            assert rhombus.contains_point(vertex)

        assert not rhombus.contains_point(Vector(rhombus.midpoint.x+3, rhombus.midpoint.y))
        assert not rhombus.contains_point(Vector(rhombus.midpoint.x-3, rhombus.midpoint.y))
        assert not rhombus.contains_point(Vector(rhombus.midpoint.x, rhombus.midpoint.y+3))
        assert not rhombus.contains_point(Vector(rhombus.midpoint.x, rhombus.midpoint.y-3))

        # a point just barely past a vertex
        assert not rhombus.contains_point(rhombus.midpoint + (vertices[0] - rhombus.midpoint)*1.01)

        midpoint = vertices[0] + (vertices[1] - vertices[0])/2
        assert rhombus.contains_point(midpoint)

        # a point just barely past the midpoint of an edge
        assert not rhombus.contains_point(rhombus.midpoint + (midpoint - rhombus.midpoint)*1.01)

        # a point just barely past the vertex
        assert not rhombus.contains_point(vertices[0] + (vertices[0] - rhombus.midpoint)*1.01)

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
            edge_midpoint = vertices[0] + (vertices[0] - vertices[1])/2

            edge_rhombus = tiling.rhombus_at_point(edge_midpoint)
            assert edge_rhombus.contains_point(edge_midpoint)

            for vertex in rhombus.vertices():
                rhombus = tiling.rhombus_at_point(vertex)
                assert rhombus.contains_point(vertex)
