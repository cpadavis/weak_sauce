"""
In here goes the various types of sources you can use to illuminate an image.

Nearly all of these should use the deposit grid method (to be incorporated in here, too) to deposit, but there is also an exact uniform illumination method as well

self.fluxes = self.source(self.vertices, self.fluxes, **kwargs)

TODO: Add method for adding sources
TODO: If there are greater than some number of vertex columns or rows, make sure you don't plot the vertex lines by default (except for plot_vertices, of course)
"""

import numpy as np
import matplotlib.pyplot as plt

from weak_sauce.shifted_cmap import shiftedColorMap
from weak_sauce.adaptive_moments.psf_evaluator import Moment_Evaluator

def zero_flux_func(vertices, **kwargs):
    return np.zeros(np.array(vertices.shape[:2]))

def init_grid(num_x=11, max_x=None, min_x=0,
              num_y=None, max_y=None, min_y=None,
              flux_func=zero_flux_func, **kwargs):
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
        # default to max_x
        max_y = max_x
    if type(min_y) == type(None):
        # default to min_x (which defaults to 0)
        min_y = min_x

    # create grid. Grid is [x,y] ordered!!!!!
    vertices = np.dstack(np.meshgrid(
        np.linspace(min_x, max_x, num_x, endpoint=True),
        np.linspace(min_y, max_y, num_y, endpoint=True),
        indexing='ij'))

    # create centroids
    centroids = vertex_centroids(vertices)

    # create fluxes
    #check if input flux_func is a function without needing galsim
    if type(flux_func) == type(np.sum):
        fluxes = flux_func(centroids, **kwargs)
    else:
        try:
            import galsim
        except ImportError:
            print 'Your flux_func is not a python function, and galsim is not installed! You should either:'
            print '1. install galsim and use galsim objects as flux_funcs'
            print '2. define your own flux_func as a python function'
            raise

        if issubclass(type(flux_func),eval('galsim.GSObject')):
            fluxes = np.zeros((centroids.shape[0],centroids.shape[1]))
            for i in np.arange(centroids.shape[0]):
                for j in np.arange(centroids.shape[1]):
                    fluxes[i,j] = flux_func.xValue(centroids[i,j,0],centroids[i,j,1])
        else:
            raise TypeError('youve got galsim installed, but your flux_func still wasnt a function or a galsim GSObject')

    return vertices, centroids, fluxes

def vertex_centroids(vertices):
    # this way we can have an array of centers of pixels of same shape as
    # the fluxes

    # super fugs
    x = vertices[:, :, 0]
    y = vertices[:, :, 1]
    x0 = x[:-1, :-1]
    x1 = x[1:, :-1]
    x2 = x[1:, 1:]
    x3 = x[:-1, 1:]
    cx = np.mean([x0, x1, x2, x3], axis=0)

    y0 = y[:-1, :-1]
    y1 = y[1:, :-1]
    y2 = y[1:, 1:]
    y3 = y[:-1, 1:]
    cy = np.mean([y0, y1, y2, y3], axis=0)

    centroids = np.dstack([cx, cy])

    return centroids

def check_vertices(vertices):

    x = vertices[:, :, 0]
    y = vertices[:, :, 1]
    x0 = x[:-1, :-1]
    x1 = x[1:, :-1]
    x2 = x[1:, 1:]
    x3 = x[:-1, 1:]

    y0 = y[:-1, :-1]
    y1 = y[1:, :-1]
    y2 = y[1:, 1:]
    y3 = y[:-1, 1:]

    # the deposit code needs clockwise convex polygons
    # this can be tested with cross product: need that cross product to be the
    # same sign and positive.
    dx01 = x1 - x0
    dx12 = x2 - x1
    dx23 = x3 - x2
    dx30 = x0 - x3
    dy01 = y1 - y0
    dy12 = y2 - y1
    dy23 = y3 - y2
    dy30 = y0 - y3

    cp = np.array([dx01 * dy12 - dx12 * dy01,
                   dx12 * dy23 - dx23 * dy12,
                   dx23 * dy30 - dx30 * dy23,
                   dx30 * dy01 - dx01 * dy30])
    if np.any(cp <= 0):
        # return that we failed
        return False
    # success means we're good!
    return True

class Source(object):
    """
    Class shell.
    """

    def __init__(self, num_x, **kwargs):
        self.vertices, self.centroids, self.fluxes = init_grid(num_x, **kwargs)

        # this makes it easier to convert
        self.r0 = self.vertices[0, 0]
        self.r1 = self.vertices[1, 1]
        self.x_min = self.vertices[:, :, 0].min()
        self.y_min = self.vertices[:, :, 1].min()
        self.x_max = self.vertices[:, :, 0].max()
        self.y_max = self.vertices[:, :, 1].max()

        self.psf_evaluator = Moment_Evaluator()

    def check_vertices(self):
        return check_vertices(self.vertices)

    def update_centroids(self):
        # if you modify the vertices, you should update the centroids, too!
        self.centroids = vertex_centroids(self.vertices)
        self.r0 = self.vertices[0, 0]
        self.r1 = self.vertices[1, 1]
        self.x_min = self.vertices[:, :, 0].min()
        self.y_min = self.vertices[:, :, 1].min()
        self.x_max = self.vertices[:, :, 0].max()
        self.y_max = self.vertices[:, :, 1].max()

    def evaluate_psf(self):
        # evaluate moments of fluxes image naievely in pixel coordinates
        return self.psf_evaluator(self.fluxes)

    def evaluate_sex(self):
        """
        Run sextractor on a numpy array (or weak_sauce source) to do photometry.
        """
        import sewpy, os
        sew = sewpy.SEW(params=["X_IMAGE", "Y_IMAGE", "FLUX_AUTO", "FLUX_ISO", "FLUX_ISOCOR", "FLAGS"],
                config={"DETECT_MINAREA":10, "PHOT_FLUXFRAC":"0.3, 0.5, 0.8"}, sexpath='/afs/slac/g/ki/software/local/bin/sex')
        from astropy.io import fits
        hdu = fits.PrimaryHDU(self.fluxes)
        hdulist = fits.HDUList([hdu])
        try:
            hdulist.writeto('iamlame.fits')
        except IOError:
            print 'you already have a file called iamlame.fits in this directory--you should rename it and then seek help...'
            raise
        out = sew('iamlame.fits')
        os.remove('iamlame.fits')
        return out["table"] # this is an astropy table.

    def plot(self, ZZ, XX=None, YY=None, fig=None, ax=None,
             pcolormesh_kwargs_in={}):

        pcolormesh_kwargs = {'linewidths': 0, 'edgecolors': 'k'}
        # if we have a small number of points, put the edges back in
        if np.sqrt(np.prod(np.shape(ZZ))) < 30:
            pcolormesh_kwargs['linewidths'] = 2
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
            # this means we have uniform flux
            cmap = plt.cm.gray_r
            # this sets the faces transparent
            vmin = np.median(ZZ) + 1
            vmax = np.median(ZZ) + 2
            cmap.set_under(alpha=0)
        elif midpoint <= 0 or midpoint >= 1:
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

    def plot_pixel_grid(self, fig=None, ax=None):
        # assume each pixel is equal size
        ZZ = self.fluxes
        fig, ax = self.plot(ZZ, fig=fig, ax=ax)
        return fig, ax

    def plot_naieve_grid(self, fig=None, ax=None):
        print('Warning! plot_naieve_grid deprecated! Use plot_pixel_grid!')
        return self.plot_pixel_grid(fig=fig, ax=ax)

    def plot_vertices(self, fig=None, ax=None):
        # just the distorted grid
        XX = self.vertices[:, :, 0]
        YY = self.vertices[:, :, 1]
        ZZ = np.zeros(np.array(XX.shape) - 1)
        pcolormesh_kwargs_in = {'linewidths': 2}
        fig, ax = self.plot(ZZ, XX, YY, fig=fig, ax=ax,
                            pcolormesh_kwargs_in=pcolormesh_kwargs_in)
        return fig, ax

"""
# example source with concave vertex
from weak_sauce.sources import Source
source = Source(num_x=4)
source.vertices[2,2] = [1.1, 1.1]
assert check_vertices(source.vertices) == False
source = Source(num_x=4)
source.vertices[2,2] = [0.8, 0.8]
vertices = source.vertices
"""
