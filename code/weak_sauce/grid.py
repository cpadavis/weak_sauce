"""
In here goes the grid class.
"""

import numpy as np


class MoveableGrid(object):
    """A mapper between two sets of grids -- say pixel and focal plane.

    Vertices is an M x N x 2 array. As long as pixels don't cross (and they
    better not!), all boxes are defined as box[m, n] having vertices v[m, n],
    v[m + 1, n], v[m, n + 1], and v[m + 1, n + 1]. The 2 are for x and y.

    Terminology: real grids are the actual coordinates of the vertices.
                 ideal grids are from the indices of the fluxes
    """

    def __init__(self, source, mover, **kwargs):
        self.source = source

        # class objects with __call__ methods
        self.mover = mover

    def update_vertices(self, **kwargs):
        # some method that modifies self.vertices
        self.source.vertices = self.mover.move_vertices(
                self.source.vertices, self.source.fluxes, **kwargs)

    def update_fluxes(self, **kwargs):
        # some method that updates fluxes. This can depend on the vertices
        # and also on some external function
        # This basically calls self.source function and adds it to the
        # fluxes
        self.source.fluxes = self.mover.move_fluxes(
                self.source.vertices, self.source.fluxes, **kwargs)

    def step(self, **kwargs):
        # self.update_fluxes(**kwargs)
        # self.update_vertices(**kwargs)
        self.source = self.mover(self.source, **kwargs)

    # wrap to the source object
    def evaluate_moments(self):
        # evaluate moments of fluxes image naievely.
        return self.source.evaluate_moments()

    def plot_real_grid(self, fig=None, ax=None):
        return self.source.plot_real_grid(fig=fig, ax=ax)

    def plot_naieve_grid(self, fig=None, ax=None):
        return self.source.plot_naieve_grid(fig=fig, ax=ax)

    def plot_vertices(self, fig=None, ax=None):
        return self.plot_vertices.plot_naieve_grid(fig=fig, ax=ax)

