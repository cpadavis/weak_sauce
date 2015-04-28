import numpy as np
import ctypes
from subprocess import check_output
"""
TODO: Loop over all pixels. Deposit correctly.
TODO: Make into pythonic class so that we can ignore all this stuff.
"""

#kludge-ily make code portable
basedir_name = check_output(["git", "rev-parse", "--show-toplevel"])
r2d_libpath = basedir_name[:-1]+'/code/r3d/r2d_lib.so' #cut newline char
r2d = ctypes.CDLL(r2d_libpath)


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


r2d.r2du_faces_from_verts.argtypes = [ctypes.POINTER(type((r2d_rvec2 * 4)())), r2d_int,
                                      ctypes.POINTER(type((r2d_plane * 4)()))]




def areas(vertices):
    # flux in a box is based on total area of box
    x, y = vertices
    # each box has (e.g. x) x[i,i] x[i,i+1], x[i+1,i], x[i+1,i+1]
    fluxes = 0.5 * np.abs((x[:-1, :-1] - x[1:, 1:]) *
                          (y[:-1, 1:] -  y[1:, :-1]) -
                          (x[:-1, 1:] -  x[1:, :-1]) *
                          (y[:-1, :-1] - y[1:, 1:]))
    return fluxes

def get_box_index_bounds(vertices):
    # assume vertices are array
    xmin = np.int(np.floor(np.min(vertices[:,0])))
    ymin = np.int(np.floor(np.min(vertices[:,1])))
    xmax = np.int(np.ceil(np.max(vertices[:,0])))
    ymax = np.int(np.ceil(np.max(vertices[:,1])))

    # TODO: add case if max is outside grid?
    return xmin, ymin, xmax, ymax


def overlap_pixel(vertices):
    # get the overlap deposits of a bounding box around the vertices
    dest_window = np.ceil(np.max(vertices, axis=0)).astype(int).tolist()

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

def deposit_source_onto_grid(source, vertices, fluxes):
    # vertices are in units of pixels on the source grid
    # some code to check that fluxes and areas are the same shape
    # and that vertices is the right shape in relation to those
    Nx, Ny = fluxes.shape
    x, y = vertices
    for i in xrange(Nx):
        for j in xrange(Ny):
            #print(i, j)
            vertices_ij = np.array([[x[i, j], y[i, j]],
                                    [x[i, j + 1], y[i, j + 1]],
                                    [x[i + 1, j + 1], y[i + 1, j + 1]],
                                    [x[i + 1, j], y[i + 1, j]],
                                    ])

            vertices_ij_floored = vertices_ij - np.floor(np.min(vertices_ij,
                                                                axis=0))
            xmin, ymin, xmax, ymax = get_box_index_bounds(vertices_ij)
            source_ij = source[ymin:ymax, xmin:xmax]
            overlap_ij = overlap_pixel(vertices_ij_floored).T
            # if j == 1:
            #     import ipdb; ipdb.set_trace()
            fluxes[i, j] += np.sum(overlap_ij * source_ij)

    return fluxes

if __name__ == '__main__':

    # pretend do devon's code?

    # clockwise
    verts_py_in = np.array([[1.5, 2], [2.5, 2], [2.5, 4.1], [2, 2.5]])
    verts_py_in = np.array([[0, 0], [3.2, 0], [2, 2.2], [0.1, 2]]) + np.array([2,3])
    # set the array to minimal pixel range
    verts_py = (verts_py_in - np.floor(np.min(verts_py_in, axis=0))).tolist()

    x, y = np.array(verts_py).T
    vert_area = 0.5 * np.abs((x[2] - x[0]) * (y[3] - y[1]) -
                             (x[3] - x[1]) * (y[2] - y[0]))
    inv_vert_area = 1. / vert_area

    dest_window = np.ceil(np.max(verts_py, axis=0)).astype(int).tolist()

    print(verts_py)
    print(dest_window)

    dest_window_c = r2d_rvec2(*dest_window)
    info = r2d.r2d_init(np.prod(dest_window))

    nverts = r2d_int(4)  # assume always 4

    """
    use dest grid for each pixel as a buffer containing the fractional pixel
    deposits.

    rasterize quad returns in each dest_grid coordinate the fraction of the area of the overlapping pixel (for rasterize_quad(1, faces, ibounds))
    """
    # need to figure out the dims for a given set of vertices
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
        verts[i] = r2d_rvec2(*verts_py[i])
    # buffers
    faces = (r2d_plane * 4)()

    # r2d.r2du_faces_from_verts(r2d_rvec2* verts, r2d_int nverts, r2d_plane* faces)
    r2d.r2du_faces_from_verts.argtypes = [ctypes.POINTER(type(verts)), r2d_int,
                                          ctypes.POINTER(type(faces))]

    r2d.r2du_faces_from_verts(ctypes.pointer(verts), nverts, ctypes.pointer(faces))

    info3 = r2d.r2d_rasterize_quad(r2d_real(inv_vert_area), faces, ibounds)

    r2d.r2d_finalize()

    # print the verts
    for i in xrange(len(verts)):
        print(i, verts[i].x, verts[i].y)

    # print the faces
    for i in xrange(len(faces)):
        # what do I expect?
        print(i, faces[i].n.x, faces[i].n.y, faces[i].d)

    # print the grid
    # for i in xrange(len(dest_grid)):
    #     print(i, dest_grid_c[i])
    for i in xrange(dest_grid.size):
        x = int(i / dest_dims[1])
        y = int(i % dest_dims[1])
        dest_grid[x, y] = dest_grid_c[i]
    print(dest_grid)




