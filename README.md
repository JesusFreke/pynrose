##### pynrose - P3 Penrose Tiling Generator

This is a python library and stand-alone program to generate P3 penrose
tilings.

### Stand-alone program

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

##### Getting started
1. Download the ptgen.jar
2. `java -jar ptgen.jar > tiling.svg` to generate a basic 10mm x 10mm tiling
3. `java -jar ptgen-jar --help` to see what other options are available.

### Library

This also serves as a library that can be used to generate P3 penrose
tilings for any other use you may have.


##### Generation Algorithm
The generation is based on the de Bruijn method, where there are 5 families
of equally spaced parallel lines, and each line intersection represents a
rhombus in the penrose tiling.

You can read more
[here](https://www.mathpages.com/home/kmath621/kmath621.htm) and
[here](http://www.ams.org/publicoutreach/feature-column/fcarc-ribbons)
