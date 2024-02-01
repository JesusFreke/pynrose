# pynrose - P3 Penrose Tiling Generator

This is a python library and stand-alone program to generate P3 penrose
tilings.

## Installation

```pip install pynrose```

## Stand-alone program

As a stand-alone program, this is able to generate P3 penrose tilings
to SVG. It supports a two modes of operation.

The SVG mode outputs each rhombus as a separate closed path, with
different styles applied to thick and thin rhombii. This is intended
for visual/display applications.

The SVGLINE mode outputs each rhombus edge as a separate path, and
ensures the edges are deduplicated, to avoid re-cutting/etching the
same line multiple times, for CNC/laser etching types of applications. 

In both modes, the tiling can be split up into multiple smaller tilings
that can be recombined with no overlaps or missing rhombii. This allows
a tiling to be split up into multiple smaller parts, based on the
constraints of the manufacturing process, and then assembled into a single
large tiling.

## CLI Usage

To generate a random penrose tiling as an SVG 

```pynrose > tiling.svg```

To explore more of the options that are available, you can use --help or -?

```pynrose --help```

## API Usage

This also exposes an API that can be used to generate P3 penrose tilings programatically
for any other use you may have.

The API documentation can be found at https://pynrose.readthedocs.io

```python
from pynrose import Tiling, Grid, Vector

# generate a new tiling with random offsets
tiling = Tiling()
grid = Grid(Vector(0, 0), Vector(20, 20))

# iterate over all rhombii whose midpoints are in the grid cell from (0, 0) to (20, 20)
for rhombus in tiling.rhombii(grid.cell(0, 0)):
    print(rhombus.vertices)
```

## Generation Algorithm
The generation algorithm is based on the de Bruijn method, where there are 5 families
of equally spaced parallel lines, and each line intersection represents a
rhombus in the penrose tiling.

You can read more at

* https://www.mathpages.com/home/kmath621/kmath621.htm
* http://www.ams.org/publicoutreach/feature-column/fcarc-ribbons