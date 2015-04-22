"""

Fundamentally this is what we have:

    - Two coordinate systems: the focal plane coordinate system and the pixel
      coordinate system.
    - The pixel coordinate system's vertices are defined in relation to the
      focal plane coordinate system. (I guess in the future it would be nice to
      be able to invert this as well.)
    - An incoming flux model defined against the focal plane coordinate system.
      That is, each 'voxel' in the focal plane coordinate system is defined by
      this model.
    - A set of values for the boxes as defined by the pixel coordinate system's
      vertices.

We can then iterate another wave of the flux model, which can then adjust the
pixel vertices. Rinse repeat.

At any given point we should be able to take those flux values and make a plot
wherein we pretend that the different boxes are all square, equally-sized
pixels.  As long as pixel vertices don't cross, this should be well-defined.


TODO: Incorporate ability to use non-uniform illumination
TODO: Figure out better way to define pixel size scale against perturbations
TODO: Incorporate timestep

"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def stationary(vertices):
    # a type of mover. Just returns the vertices
    return vertices


def uniform_flux(vertices):
    # flux in a box is based on total area of box
    x, y = vertices
    # each box has (e.g. x) x[i,i] x[i,i+1], x[i+1,i], x[i+1,i+1]
    fluxes = 0.5 * np.abs((x[:-1, :-1] - x[1:, 1:]) *
                          (y[:-1, 1:] -  y[1:, :-1]) -
                          (x[:-1, 1:] -  x[1:, :-1]) *
                          (y[:-1, :-1] - y[1:, 1:]))
    return fluxes

def treering(vertices, center, wavelength, amplitude=1, phi=0):
    x, y = vertices
    x0, y0 = center

    r = np.sqrt((x - x0) ** 2 + (y - y0) ** 2)
    theta = np.arctan2((y - y0), (x - x0))

    # displace along radius with some sort of wavelength
    r += amplitude * np.sin(2 * np.pi * r / wavelength + phi)
    x = r * np.cos(theta) + x0
    y = r * np.sin(theta) + y0

    return [x, y]


class MoveableGrid(object):
    """A mapper between two sets of grids -- say pixel and focal plane.

    Vertices is an M x N array. As long as pixels don't cross (and they better
    not!), all boxes are defined as box[m, n] having vertices v[m, n], v[m + 1,
    n], v[m, n + 1], and v[m + 1, n + 1].

    Terminology: real grids are the actual coordinates of the vertices.
                 ideal grids are from the indices of the fluxes
    """

    def __init__(self, vertices, fluxes, source=uniform_flux, mover=stationary):
        # TODO: make sure vertices and fluxes are the right shapes
        self.vertices = vertices
        self.fluxes = fluxes
        self.source = source
        self.mover = mover

    def move_vertices(self, **kwargs):
        # some method that modifies self.vertices
        self.vertices = self.mover(self.vertices, **kwargs)

    def update_fluxes(self, **kwargs):
        # some method that updates fluxes. This can depend on the vertices
        # and also on some external function
        # This basically calls self.illumination function and adds it to the
        # fluxes
        self.fluxes += self.source(self.vertices, **kwargs)

    def step(self, **kwargs):
        self.update_fluxes(**kwargs)
        self.move_vertices(**kwargs)

    def plot_grid(self, ZZ, XX=None, YY=None, fig=None, ax=None,
            pcolormesh_kwargs_in={}):
        # TODO: add kwargs argument for pcolormesh
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
            vmin = np.median(ZZ) + 1
            vmax = np.median(ZZ) + 2
            cmap.set_under(alpha=0)
        elif midpoint <= 0:
            cmap = plt.cm.RdBu_r#plt.cm.Reds
        else:
            cmap = shiftedColorMap(plt.cm.RdBu_r, midpoint=midpoint)
        if type(XX) == type(None):
            IM = ax.pcolormesh(ZZ, cmap=cmap, vmin=vmin, vmax=vmax, **pcolormesh_kwargs)
        else:
            IM = ax.pcolormesh(XX, YY, ZZ, cmap=cmap, vmin=vmin, vmax=vmax, **pcolormesh_kwargs)
        if not np.all(ZZ == np.median(ZZ)):
            fig.colorbar(IM, ax=ax)
        else:
            #import ipdb; ipdb.set_trace()
            #IM.set_alpha(0.1)
            #IM.set_facecolor((0,1,0,0))
            pass
        if type(XX) != type(None):
            ax.set_xlim(XX.min(), XX.max())
        else:
            ax.set_xlim(0, ZZ.shape[1])
        if type(YY) != type(None):
            ax.set_ylim(YY.min(), YY.max())
        else:
            ax.set_ylim(0, ZZ.shape[0])
        return fig, ax

    def plot_real_grid(self, fig=None, ax=None):
        # use vertices
        XX, YY = self.vertices
        ZZ = self.fluxes
        fig, ax = self.plot_grid(ZZ, XX, YY, fig=fig, ax=ax)
        return fig, ax

    def plot_naieve_grid(self, fig=None, ax=None):
        # assume each pixel is equal size
        ZZ = self.fluxes
        num_x, num_y = ZZ.shape
        num_x += 1
        num_y += 1
        x, y = self.vertices
        XX, YY = np.meshgrid(np.linspace(x.min(), x.max(), num_x, endpoint=True),
                             np.linspace(y.min(), y.max(), num_y, endpoint=True))
        fig, ax = self.plot_grid(ZZ, fig=fig, ax=ax)#, XX, YY, fig=fig, ax=ax)
        return fig, ax

    def plot_vertices(self, fig=None, ax=None):
        # just the distorted grid
        XX, YY = self.vertices
        ZZ = np.zeros(np.array(XX.shape) - 1)
        fig, ax = self.plot_grid(ZZ, XX, YY, fig=fig, ax=ax)
        return fig, ax

    def plot_source(self, fig=None, ax=None):
        # assume each pixel is equal size
        ZZ = self.source
        num_x, num_y = ZZ.shape
        num_x += 1
        num_y += 1
        XX, YY = np.meshgrid(np.linspace(0, num_x - 1, num_x, endpoint=True),
                             np.linspace(0, num_y - 1, num_y, endpoint=True))
        fig, ax = self.plot_grid(ZZ, fig=fig, ax=ax, pcolormesh_kwargs_in={'linewidths': 0})
        return fig, ax


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
Taken from

https://github.com/olgabot/prettyplotlib/blob/master/prettyplotlib/colors.py

which makes beautiful plots by the way


    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

if __name__ == '__main__':

    # demonstrate small perturbations to pixel
    num = 21
    x, y = np.meshgrid(np.linspace(-1, 1, num, endpoint=True),
                       np.linspace(-1, 1, num, endpoint=True))
    flux_norm = (x[0,0] - x[0,1]) ** 2
    # move the vertices around
    y[10:] *= 2
    distortion = 0.2 * 2. / num
    x += np.random.random(x.shape) * 2 * distortion - distortion
    y += np.random.random(y.shape) * 2 * distortion - distortion

    fluxes = (uniform_flux([x, y]) - flux_norm) / flux_norm

    grid = MoveableGrid([x, y], fluxes)
    fig, ax = grid.plot_vertices()
    fig, ax = grid.plot_real_grid()
    fig, ax = grid.plot_naieve_grid()


    # demonstrate tree rings
    num = 21
    x, y = np.meshgrid(np.linspace(-1, 1, num, endpoint=True),
                       np.linspace(-1, 1, num, endpoint=True))
    L = (x[0,0] - x[0,1])
    flux_norm = L ** 2
    # move the vertices around
    x, y = treering([x, y], [0.5, 0.5], 10 * L, 0.5 * L, np.pi)
    fluxes = (uniform_flux([x, y]) - flux_norm) / flux_norm

    grid = MoveableGrid([x, y], fluxes)
    fig, ax = grid.plot_vertices()
    fig, ax = grid.plot_real_grid()
    fig, ax = grid.plot_naieve_grid()


    # try to do devon's fig 8 just using the exact integration
    num = 32
    x, y = np.meshgrid(np.linspace(-0.5, 0.5, num, endpoint=True),
                       np.linspace(-0.5, 0.5, num, endpoint=True))
    L = (x[0,0] - x[0,1])
    flux_norm = L ** 2

    r = np.sqrt(x ** 2 + y ** 2)
    theta = np.arctan2(y, x)
    conds = (r / 0.45) < 1
    # sign difference? or just I'm finding area and he's finding density
    theta[conds] = theta[conds] + np.pi / 2 * ((r[conds] / 0.45) ** 2 - 1) ** 2
    r[conds] = r[conds] * (1 + 0.2 * (r[conds] / 0.45 - 1) ** 2)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    fluxes = (uniform_flux([x, y])) / flux_norm
    grid = MoveableGrid([x, y], fluxes)
    fig, ax = grid.plot_vertices()
    fig, ax = grid.plot_real_grid()
    fig, ax = grid.plot_naieve_grid()


    from test_r2d import deposit_source_onto_grid
    # TODO: get vertices in focal plane coords, not pixel coords. Deal with the
    # transformation.
    # TODO: Fix this x, y ordering of stuff. Meshgrids produce y-ordered
    # indexing, while the code assumes x ordering.
    xnum = 20
    ynum = 25
    xmax = 50
    ymax = 70
    # x, y = np.meshgrid(np.linspace(0, xmax, xnum, endpoint=True),
    #                    np.linspace(0, np.sqrt(ymax), ynum, endpoint=True) ** 2)
    x, y = np.meshgrid(np.linspace(0, xmax, xnum, endpoint=True),
                       np.linspace(0, ymax, ynum, endpoint=True))
    # rotate the pixels
    r = np.sqrt((x - 0.5 * xmax) ** 2 + (y - 0.5 * ymax) ** 2)
    theta = np.arctan2(y - 0.5 * ymax, x - 0.5 * xmax)
    theta -= np.pi / 4
    # twist them
    theta += np.pi / 4 * ( (r / r.max()) ** 2)
    # stretch them
    r = r ** 2
    x = r * np.cos(theta) + 0.5 * xmax
    x -= x.min()
    x *= xmax / x.max()
    y = r * np.sin(theta) + 0.5 * ymax
    y -= y.min()
    y *= ymax / y.max()
    vertices = [x, y]
    areas = uniform_flux(vertices)
    ## Uniform Illumination.
    source = np.ones([ymax + 1, xmax + 1])
    # xx, yy = np.meshgrid(np.arange(0, xmax + 1),
    #                      np.arange(0, ymax + 1))
    # source = (xx - xmax * 0.25) ** 2 + 5 * (yy - ymax * 0.25) ** 2 + -4 * (xx - xmax * 0.3) * (yy - ymax * 0.75)
    # source /= source.sum()
    # source -= source.mean()

    fluxes_zeros = np.zeros(areas.shape)
    fluxes_prime = deposit_source_onto_grid(source, vertices, fluxes_zeros)
    grid = MoveableGrid([x, y], fluxes_prime, source=source)
    fig, ax = grid.plot_vertices()
    fig, ax = grid.plot_real_grid()
    fig, ax = grid.plot_naieve_grid()

    # confirm with uniform illumination
    grid = MoveableGrid([x, y], uniform_flux([x, y]))
    fig, ax = grid.plot_vertices()
    fig, ax = grid.plot_real_grid()
    fig, ax = grid.plot_naieve_grid()

    # Cray Illumination Pattern
    xx, yy = np.meshgrid(np.arange(0, xmax + 1),
                         np.arange(0, ymax + 1))
    source = (xx - xmax * 0.25) ** 2 + 5 * (yy - ymax * 0.25) ** 2 + -4 * (xx - xmax * 0.3) * (yy - ymax * 0.75)
    source /= source.sum()
    source -= source.mean()

    fluxes_zeros = np.zeros(areas.shape)
    fluxes_prime = deposit_source_onto_grid(source, vertices, fluxes_zeros)
    grid = MoveableGrid([x, y], fluxes_prime, source=source)
    fig, ax = grid.plot_source()
    fig, ax = grid.plot_vertices(fig=fig, ax=ax)
    fig, ax = grid.plot_real_grid()
    fig, ax = grid.plot_naieve_grid()

    plt.show()

    #import ipdb; ipdb.set_trace()
