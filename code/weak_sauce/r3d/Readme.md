# r3d

By Devon Powell

## About:

Functions for exactly voxelizing 3D primitives (tetrahedra) onto a 3D Cartesian grid including
higher-order moments, as described in [Powell & Abel (2014)](http://arxiv.org/abs/1412.4941).

Also includes utilities for carrying out basic clipping operations on convex polyhedra in a fast,
geometrically robust way. This software thus forms the basis for an exact general remeshing scheme.

This functionality is provided for 2D as well.

## Features:

- Voxelizes a tetrahedron onto a regular grid by calculating the exact coordinate moments (1, x, y,
  z, x^2, y^2, z^2, xy, xz, yz) of the intersection between the tetrahedron and each voxel.
- All declarations and documentation are located in r3d.h and r2d.h.

## TODO:

- Add generalized functions for voxelizing arbitrary polyhedra, or even polygons embedded in R^3.

## Usage:

- To build, type

`make`

- To compile into your code,

`#include <r3d.h>`

- To link,

`-lr3d`

- - - 

Copyright (C) 2014 Stanford University. See License.txt for more information.

