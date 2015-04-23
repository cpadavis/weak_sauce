"""
In here goes the various types of sources you can use to illuminate an image.

Nearly all of these should use the deposit grid method (to be incorporated in here, too) to deposit, but there is also an exact uniform illumination method as well

self.fluxes = self.source(self.vertices, self.fluxes, **kwargs)


"""

import numpy as np
import matplotlib.pyplot as plt

from weak_sauce.shifted_cmap import shiftedColorMap
from weak_sauce.adaptive_moments.psf_evaluator import Moment_Evaluator


def init_grid(num_x=11, max_x=None, min_x=0,
              num_y=None, max_y=None, min_y=None):
    """
    z[i, j] returns x_ij, y_ij
    f[i, j] returns f_ij
    """
    # set out all the parameters
    if type(max_x) == type(None):
        # default to pixel coordinates
        max_x = num_x - 1
    if type(num_y) == type(None):
        # default to square grid
        num_y = num_x
    if type(max_y) == type(None):
        # default to pixel coordinates
        max_y = num_y - 1
    if type(min_y) == type(None):
        # default to min_x (which defaults to 0)
        min_y = min_x

    # create grid. Grid is [x,y] ordered!!!!!
    vertices = np.dstack(np.meshgrid(
        np.linspace(min_x, max_x, num_x, endpoint=True),
        np.linspace(min_y, max_y, num_y, endpoint=True),
        indexing='ij'))

    # create fluxes
    fluxes = np.zeros(np.array(vertices.shape[:2]) - 1)

    return vertices, fluxes

def vertex_centroids(vertices):
    # this way we can have an array of centers of pixels of same shape as
    # the fluxes

    # super fugs
    x = vertices[:, :, 0]
    y = vertices[:, :, 1]
    x0 = x[:-1, :-1]
    x1 = x[:-1, 1:]
    x2 = x[1:, 1:]
    x3 = x[1:, :-1]
    cx = np.mean([x0, x1, x2, x3], axis=0)

    y0 = y[:-1, :-1]
    y1 = y[:-1, 1:]
    y2 = y[1:, 1:]
    y3 = y[1:, :-1]
    cy = np.mean([y0, y1, y2, y3], axis=0)

    centroids = np.dstack([cx, cy])

    return centroids

class Source(object):
    """
    Class shell.
    """

    def __init__(self, init_grid_func=init_grid, **kwargs):
        self.vertices, self.fluxes = init_grid_func(**kwargs)
        self.centroids = vertex_centroids(self.vertices)

        # this makes it easier to convert
        self.r0 = self.vertices[0, 0]
        self.r1 = self.vertices[1, 1]

        self.moment_evaluator = Moment_Evaluator()

    def update_centroids(self):
        # if you modify the vertices, you should update the centroids, too!
        self.centroids = vertex_centroids(self.vertices)


    def evaluate_moments(self):
        # evaluate moments of fluxes image naievely in pixel coordinates
        return self.moment_evaluator(self.fluxes)

    def plot(self, ZZ, XX=None, YY=None, fig=None, ax=None,
             pcolormesh_kwargs_in={}):

        pcolormesh_kwargs = {'linewidths': 2, 'edgecolors': 'k'}
        pcolormesh_kwargs.update(pcolormesh_kwargs_in)
        if type(fig) == type(None):
            fig, ax = plt.subplots()
        # figure out shifting the colormap
        b = np.max(ZZ)
        a = np.min(ZZ)
        c = 0
        midpoint = (c - a) / (b - a)
        vmin = a
        vmax = b
        if np.all(ZZ == np.median(ZZ)):
            cmap = plt.cm.gray_r
            # this sets the faces transparent
            vmin = np.median(ZZ) + 1
            vmax = np.median(ZZ) + 2
            cmap.set_under(alpha=0)
        elif midpoint <= 0:
            cmap = plt.cm.RdBu_r#plt.cm.Reds
        else:
            cmap = shiftedColorMap(plt.cm.RdBu_r, midpoint=midpoint)

        if type(XX) == type(None):
            # make an XX and YY grid
            num_x, num_y = ZZ.shape
            XX, YY = np.meshgrid(np.linspace(0, num_x, num_x + 1, endpoint=True),
                                 np.linspace(0, num_y, num_y + 1, endpoint=True),
                                 indexing='ij')
        IM = ax.pcolormesh(XX, YY, ZZ, cmap=cmap, vmin=vmin, vmax=vmax,
                           **pcolormesh_kwargs)
        if not np.all(ZZ == np.median(ZZ)):
            # don't plot a colorbar if we don't have any values!
            fig.colorbar(IM, ax=ax)

        ax.set_xlim(XX.min(), XX.max())
        ax.set_ylim(YY.min(), YY.max())

        return fig, ax

    def plot_real_grid(self, fig=None, ax=None):
        # use vertices
        XX = self.vertices[:, :, 0]
        YY = self.vertices[:, :, 1]
        ZZ = self.fluxes
        fig, ax = self.plot(ZZ, XX, YY, fig=fig, ax=ax)
        return fig, ax

    def plot_naieve_grid(self, fig=None, ax=None):
        # assume each pixel is equal size
        ZZ = self.fluxes
        fig, ax = self.plot(ZZ, fig=fig, ax=ax)
        return fig, ax

    def plot_vertices(self, fig=None, ax=None):
        # just the distorted grid
        XX = self.vertices[:, :, 0]
        YY = self.vertices[:, :, 1]
        ZZ = np.zeros(np.array(XX.shape) - 1)
        fig, ax = self.plot(ZZ, XX, YY, fig=fig, ax=ax)
        return fig, ax


class StationarySource(Source):
    """
    Idea is here you give a fixed array type object and then deposit it onto
    a grid.
    """
    def __init__(self, source_function,
                 num_x=11, max_x=None, min_x=0,
                 num_y=None, max_y=None, min_y=None):
        super(StationarySource, self).__init__(init_grid_func=init_grid,
                num_x=num_x, max_x=max_x, min_x=min_x,
                num_y=num_y, max_y=max_y, min_y=min_y)
        # some kind of way of specifying any pixel coordinate transformations?

        # now create the source array
        self.source_function = source_function
        self.fluxes = self.source_function(self.centroids)

