Usage
=====

.. _installation:

Installation
------------

.. code-block:: console

   pip install pynrose

Running from the CLI
--------------------

To generate a penrose tiling as an SVG file from the CLI, you can run

.. code-block:: console

   pynrose > penrose.svg

To explore more of the options that are available, you can use --help or -?

.. code-block:: console

   pynrose --help


Using the API
------------------

Example usage

.. code-block:: python

  from pynrose import Tiling, Grid, Vector

  # generate a new tiling with random offsets
  tiling = Tiling()
  grid = Grid(Vector(0, 0), Vector(20, 20))

  # iterate over all rhombii whose midpoints are in the grid cell from (0, 0) to (20, 20)
  for rhombus in tiling.rhombii(grid.cell(0, 0)):
    print(rhombus.vertices)
