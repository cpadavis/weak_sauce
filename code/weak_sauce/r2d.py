"""
Interfaces with the r2d C library
"""

import numpy as np
import ctypes
import os

# load library
r2d = ctypes.CDLL(os.path.dirname(os.path.realpath(__file__)) +
                  '/r3d/r2d_lib.so')

# define some types
r2d_int = ctypes.c_int
r2d_long = ctypes.c_int
#r2d_real = ctypes.c_float
r2d_real = ctypes.c_double

# define some structures

class r2d_rvec2(ctypes.Structure):
    _fields_ = [('x', r2d_real),
                ('y', r2d_real)]

class r2d_dvec2(ctypes.Structure):
    _fields_ = [('x', r2d_int),
                ('y', r2d_int)]

class r2d_plane(ctypes.Structure):
    _fields_ = [('n', r2d_rvec2),
                ('d', r2d_real)]

class r2d_vertex(ctypes.Structure):
    _fields_ = [('pos', r2d_rvec2),
                ('pnbrs', ctypes.c_ubyte * 2),
                ('fdist', r2d_real * 2),
                ('fflags', ctypes.c_ubyte)]

class r2d_info(ctypes.Structure):
    _fields_ = [('good', r2d_int),
                ('errmsg', ctypes.POINTER(ctypes.c_char)),
                ('vtot', r2d_real),
                ('vox_min', r2d_real),
                ('vox_max', r2d_real)]

# mess with functions

# r2d.r2d_init
r2d.r2d_init.restype = r2d_info
r2d.r2d_init.argtypes = [r2d_int]


# r2d.r2d_set_dest_grid
r2d.r2d_set_dest_grid.restype = r2d_info
#r2d.r2d_set_dest_grid.argtypes = [ctypes.POINTER(r2d_real),
#                                  r2d_dvec2, r2d_rvec2]
# r2d.r2d_finalize
#r2d.r2d_finalize.restype = ctypes.c_void
#r2d.r2d_finalize.argtypes
# r2d.r2d_rasterize_quad
r2d.r2d_rasterize_quad.restype = r2d_info
r2d.r2d_rasterize_quad.argtypes = [r2d_real, (r2d_plane * 4), (r2d_dvec2 * 2)]


def get_box_index_bounds(vertices):
    # assume vertices are array
    xmin = np.int(np.floor(np.min(vertices[:,0])))
    ymin = np.int(np.floor(np.min(vertices[:,1])))
    xmax = np.int(np.floor(np.max(vertices[:,0])))
    ymax = np.int(np.floor(np.max(vertices[:,1])))

    # TODO: add case if max is outside grid?
    return xmin, ymin, xmax, ymax


def overlap_pixel(vertices, dest_window):
    # get the overlap deposits of a bounding box around the vertices
    #dest_window = np.ceil(np.max(vertices, axis=0)).astype(int).tolist()

    dest_window_c = r2d_rvec2(*dest_window)
    info = r2d.r2d_init(np.prod(dest_window))
    nverts = r2d_int(4)  # assume always 4

    dest_dims = dest_window
    dest_grid = np.zeros(np.prod(dest_dims))
    ibounds = (r2d_dvec2 * 2)()
    ibounds[0].x = 0
    ibounds[0].y = 0
    ibounds[1].x = dest_dims[0]
    ibounds[1].y = dest_dims[1]

    dest_dims_c = r2d_dvec2(*dest_dims)
    dest_grid_c = (r2d_real * len(dest_grid))(*dest_grid)
    dest_grid = dest_grid.reshape(dest_dims)
    r2d.r2d_set_dest_grid.argtypes = [ctypes.POINTER(type(dest_grid_c)),
                                      r2d_dvec2, r2d_rvec2]
    info2 = r2d.r2d_set_dest_grid(ctypes.pointer(dest_grid_c), dest_dims_c,
                                  dest_window_c)

    # put in some verts
    # verts = (r2d_vertex * 4)()
    # despite their name, they are not of the vertex class.
    verts = (r2d_rvec2 * 4)()
    for i in xrange(len(verts)):
        verts[i] = r2d_rvec2(*vertices[i])
    # buffers
    faces = (r2d_plane * 4)()

    r2d.r2du_faces_from_verts.argtypes = [ctypes.POINTER(type(verts)), r2d_int,
                                          ctypes.POINTER(type(faces))]
    r2d.r2du_faces_from_verts(ctypes.pointer(verts), nverts,
                              ctypes.pointer(faces))

    info3 = r2d.r2d_rasterize_quad(r2d_real(1), faces, ibounds)

    r2d.r2d_finalize()

    # deposit the results
    for i in xrange(dest_grid.size):
        x = int(i / dest_dims[1])
        y = int(i % dest_dims[1])
        dest_grid[x, y] = dest_grid_c[i]

    return dest_grid

def deposit(source_array, vertices):
    # vertices are in units of pixels on the source grid
    x = vertices[:, :, 0]
    y = vertices[:, :, 1]
    Nx, Ny = np.array(vertices.shape[:2]) - 1
    fluxes = np.zeros((Nx, Ny))

    # check that the vertices are within the source_array limits
    Sx, Sy = np.array(source_array.shape) + 1

    for i in xrange(Nx):
        for j in xrange(Ny):

            vertices_ij = np.array([[x[i, j], y[i, j]],
                                    [x[i + 1, j], y[i + 1, j]],
                                    [x[i + 1, j + 1], y[i + 1, j + 1]],
                                    [x[i, j + 1], y[i, j + 1]],
                                    ])

            vertices_ij_floored = vertices_ij - np.floor(np.min(vertices_ij,
                                                                axis=0))
            xmin, ymin, xmax, ymax = get_box_index_bounds(vertices_ij)
            if xmin < 0:
                xmin = 0
            if ymin < 0:
                ymin = 0
            if xmax >= source_array.shape[0]:
                xmax = source_array.shape[0] - 1
            if ymax >= source_array.shape[1]:
                ymax = source_array.shape[1] - 1
            source_ij = source_array[xmin:xmax + 1, ymin:ymax + 1]
            dest_window = [xmax-xmin + 1, ymax-ymin + 1]
            overlap_ij = overlap_pixel(vertices_ij_floored, dest_window)
            # if j == 6:
            #     # import ipdb; ipdb.set_trace()
            if np.any(np.array(overlap_ij.shape) != np.array(source_ij.shape)):
                print(i, j)
                print(overlap_ij)
                print(source_ij)
                print(vertices_ij_floored)
                print(vertices_ij)
                print(xmin, ymin, xmax, ymax)
                print(vertices.shape, source_array.shape)
                print(overlap_ij.shape, source_ij.shape)
            fluxes[i, j] += np.sum(overlap_ij * source_ij)

    return fluxes

