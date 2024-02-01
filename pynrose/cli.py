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
import argparse
import math
import random

from pynrose import *


def generate_svg(
        tiling: Tiling,
        grid: Grid,
        grid_count_x: int, grid_count_y: int,
        grid_spacing_x: float, grid_spacing_y: float,
        add_grid_box: bool):

    # the maximum amount a tile can stick out of the bounding box. Basically, half of a thin rhombus
    max_protrusion = math.sin(math.radians(72))

    view_size = Vector(
        grid.grid_size.x * grid_count_x + ((grid_count_x - 1) * grid_spacing_x) + max_protrusion * 2,
        grid.grid_size.y * grid_count_y + ((grid_count_y - 1) * grid_spacing_y) + max_protrusion * 2)

    print('<svg width="%fmm" height="%fmm" viewBox="%f %f %f %f">' % (
        view_size.x,
        view_size.y,
        grid.origin.x - max_protrusion,
        grid.origin.y - max_protrusion,
        view_size.x,
        view_size.y))
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

    for grid_x in range(grid_count_x):
        for grid_y in range(grid_count_y):
            offset = Vector(grid_x * grid_spacing_x, grid_y * grid_spacing_y)

            for rhombus in tiling.rhombii(grid.cell(grid_x, grid_y)):
                string = '<path'

                if rhombus.type() == RhombusType.THIN:
                    string += ' class="thinRhombus"'
                else:
                    string += ' class="thickRhombus"'

                string += ' id="Rhombus (%d, %d) (%d, %d) %s"' % (
                    rhombus.strip1.family.pentangle.pentangle,
                    rhombus.strip1.multiple,
                    rhombus.strip2.family.pentangle.pentangle,
                    rhombus.strip2.multiple,
                    [*rhombus.lattice_coords])

                string += ' d="M'

                for vertex in rhombus.vertices():
                    coordinate = vertex.coordinate + offset
                    string += ' %f,%f' % (coordinate.x, coordinate.y)
                string += ' z"/>'
                print(string)

            if add_grid_box:
                print('<rect x="%f" y="%f" width="%f" height="%f" class="boundingBox"/>' % (
                    grid.origin.x + grid_x * grid.grid_size.x + grid_spacing_x * grid_x,
                    grid.origin.y + grid_y * grid.grid_size.y + grid_spacing_y * grid_y,
                    grid.grid_size.x,
                    grid.grid_size.y))
    print('</svg>')


def generate_svgline(
        tiling: Tiling,
        grid: Grid,
        grid_count_x: int, grid_count_y: int,
        grid_spacing_x: float, grid_spacing_y: float,
        add_grid_box: bool):

    # the maximum amount a tile can stick out of the bounding box. Basically, half of a thin rhombus
    max_protrusion = math.sin(math.radians(72))

    view_size = Vector(
        grid.grid_size.x * grid_count_x + ((grid_count_x - 1) * grid_spacing_x) + max_protrusion * 2,
        grid.grid_size.y * grid_count_y + ((grid_count_y - 1) * grid_spacing_y) + max_protrusion * 2)

    print('<svg width="%fmm" height="%fmm" viewBox="%f %f %f %f">' % (
        view_size.x,
        view_size.y,
        grid.origin.x - max_protrusion,
        grid.origin.y - max_protrusion,
        view_size.x,
        view_size.y))
    print('<style><![CDATA[')
    print('rect.boundingBox {')
    print('    stroke: blue;')
    print('    stroke-width: .05;')
    print('    fill-opacity: 0;')
    print('    stroke-opacity: .5;')
    print('}')
    print('path.rhombusEdge {')
    print('    stroke: #000000;')
    print('    stroke-width: .01;')
    print('}')
    print(']]></style>')

    for grid_x in range(grid_count_x):
        for grid_y in range(grid_count_y):
            offset = Vector(grid_x * grid_spacing_x, grid_y * grid_spacing_y)

            for vertex1, vertex2 in tiling.rhombus_edges(grid.cell(grid_x, grid_y)):
                string = '<path class="rhombusEdge"'

                string += ' d="M'

                coordinate = vertex1.coordinate + offset
                string += ' %f,%f' % (coordinate.x, coordinate.y)

                coordinate = vertex2.coordinate + offset
                string += ' %f,%f' % (coordinate.x, coordinate.y)

                string += ' z"/>'
                print(string)

            if add_grid_box:
                print('<rect x="%f" y="%f" width="%f" height="%f" class="boundingBox"/>' % (
                    grid.origin.x + grid_x * grid.grid_size.x + grid_spacing_x * grid_x,
                    grid.origin.y + grid_y * grid.grid_size.y + grid_spacing_y * grid_y,
                    grid.grid_size.x,
                    grid.grid_size.y))
    print('</svg>')


def main():
    parser = argparse.ArgumentParser(
        description="Generates an SVG containing a penrose P3 tiling.",
        add_help=False)

    parser.add_argument(
        "-?", "--help",
        action='help',
        help='Show this help message and exit.')

    output_group = parser.add_argument_group("Output").add_mutually_exclusive_group()
    output_group.add_argument("--svg", action='store_true',
                              help="(Default)Generate a dual-color SVG with each rhombus as a separate path.")
    output_group.add_argument("--svgline", action='store_true',
                              help="Generate an SVG containing the individual, deduplicated rhombus edges.")

    grid_group = parser.add_argument_group("Grid")
    grid_group.add_argument("--minX", "-x", type=float, default=0.0, help="The minimum x value of the bounding grid.")
    grid_group.add_argument("--minY", "-y", type=float, default=0.0, help="The minimum y value of the bounding grid.")
    grid_group.add_argument("--width", "-w", type=float, default=20.0, help="The width of each grid cell.")
    grid_group.add_argument("--height", "-h", type=float, default=20.0, help="The height of each grid cell.")
    grid_group.add_argument("--count_x", "-cx", type=int, default=1,
                            help="The number of grids to generate in the x axis.")
    grid_group.add_argument("--count_y", "-cy", type=int, default=1,
                            help="The number of grids to generate in the y axis.")
    grid_group.add_argument("--grid_spacing_x", "-gx", type=float, default=2,
                            help="How much space to leave between each grid, in the x axis.")
    grid_group.add_argument("--grid_spacing_y", "-gy", type=float, default=2,
                            help="How much space to leave between each grid, in the y axis.")
    grid_group.add_argument("--add_grid_box", "-gb", action='store_true',
                            help="Add the grid bounding boxes to the svg.")

    parser.add_argument("--seed", "-s", type=int, help="The random seed to use to generate the tiling.")

    args = parser.parse_args()

    tiling = Tiling(rnd=random.Random(args.seed))
    grid = Grid(Vector(args.minX, args.minY), Vector(args.width, args.height))

    if args.svgline:
        generate_svgline(
            tiling, grid, args.count_x, args.count_y, args.grid_spacing_x, args.grid_spacing_y, args.add_grid_box)
    else:
        generate_svg(
            tiling, grid, args.count_x, args.count_y, args.grid_spacing_x, args.grid_spacing_y, args.add_grid_box)


if __name__ == "__main__":
    main()
