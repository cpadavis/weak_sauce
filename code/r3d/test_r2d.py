import numpy as np
import ctypes

"""
TODO: Loop over all pixels. Deposit correctly.
TODO: Make into pythonic class so that we can ignore all this stuff.
"""
r2d = ctypes.CDLL('/Users/cpd/Projects/weak_pads/code/r3d/r2d_lib.so')

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
# not sure about this r2dplane [4] and dvec2 ibouns[2]?
r2d.r2d_rasterize_quad.argtypes = [r2d_real, (r2d_plane * 4), (r2d_dvec2 * 2)]

# r2d.r2du_faces_from_verts(r2d_rvec2* verts, r2d_int nverts, r2d_plane* faces)
r2d.r2du_faces_from_verts.argtypes = [ctypes.POINTER(r2d_rvec2), r2d_int, ctypes.POINTER(r2d_plane)]


# pretend do devon's code?

dest_window = [3,3]
dest_window_c = r2d_rvec2(*dest_window)
info = r2d.r2d_init(np.prod(dest_window))

nverts = r2d_int(4)  # assume always 4

"""
use dest grid for each pixel as a buffer containing the fractional pixel
deposits.

rasterize quad returns in each dest_grid coordinate the fraction of the area of the overlapping pixel (for rasterize_quad(1, faces, ibounds))
"""
# need to figure out the dims for a given set of vertices
dest_dims = [3,3]
dest_grid = np.zeros(np.prod(dest_dims))
ibounds = (r2d_dvec2 * 2)()
ibounds[0].x = 0
ibounds[0].y = 0
ibounds[1].x = dest_dims[0]
ibounds[1].y = dest_dims[1]

dest_dims_c = r2d_dvec2(*dest_dims)
dest_grid_c = (r2d_real * len(dest_grid))(*dest_grid)
r2d.r2d_set_dest_grid.argtypes = [ctypes.POINTER(type(dest_grid_c)),
                                  r2d_dvec2, r2d_rvec2]
info2 = r2d.r2d_set_dest_grid(ctypes.pointer(dest_grid_c), dest_dims_c,
                              dest_window_c)


# put in some verts
# verts = (r2d_vertex * 4)()
# despite their name, they are not of the vertex class.
verts = (r2d_rvec2 * 4)()
verts_py = [[0, 0], [1.5, 0], [1.5, 1.5], [0, 1.5]]
for i in xrange(len(verts)):
    verts[i] = r2d_rvec2(*verts_py[i])
# buffers
faces = (r2d_plane * 4)()

# r2d.r2du_faces_from_verts(r2d_rvec2* verts, r2d_int nverts, r2d_plane* faces)
r2d.r2du_faces_from_verts.argtypes = [ctypes.POINTER(type(verts)), r2d_int,
                                      ctypes.POINTER(type(faces))]

r2d.r2du_faces_from_verts(ctypes.pointer(verts), nverts, ctypes.pointer(faces))

info3 = r2d.r2d_rasterize_quad(r2d_real(1.), faces, ibounds)

r2d.r2d_finalize()

# print the verts
for i in xrange(len(verts)):
    print(i, verts[i].x, verts[i].y)

# print the faces
for i in xrange(len(faces)):
    # what do I expect?
    print(i, faces[i].n.x, faces[i].n.y, faces[i].d)

# print the grid
for i in xrange(len(dest_grid)):
    print(i, dest_grid_c[i])


