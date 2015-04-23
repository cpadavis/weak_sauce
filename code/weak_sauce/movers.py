"""
Movers update flux AND vertices

TODO: Add method for adding movers
"""

import numpy as np

from weak_sauce.r2d import deposit

class Mover(object):
    """
    Class shell for movers. Currently super super simple. Just needs to have
    some kind of moving function that it calls.
    """

    def move_vertices(self, vertices, fluxes, **kwargs):
        raise NotImplementedError

    def move_fluxes(self, vertices, fluxes, **kwargs):
        raise NotImplementedError

    def move(self, source, **kwargs):
        vertices = self.move_vertices(source.vertices, source.fluxes,
                                             **kwargs)
        source.vertices = vertices
        source.update_centroids()
        fluxes = self.move_fluxes(source.vertices, source.fluxes,
                                     **kwargs)
        source.fluxes += fluxes

    def __call__(self, source, **kwargs):
        self.move(source, **kwargs)

class StationaryMover(Mover):
    """
    Don't move anything!
    """

    def move_vertices(self, vertices, fluxes, **kwargs):
        return vertices

    def move_fluxes(self, vertices, fluxes, **kwargs):
        return fluxes

class FixedIlluminationMover(StationaryMover):
    """
    Deposit some fixed but funny-shaped OTHER source
    """
    def __init__(self, stationary_source):
        self.stationary_source = stationary_source
        # self.centroids = self.stationary_source.centroids
        self.fluxes = self.stationary_source.fluxes
        self.r0 = self.stationary_source.r0
        self.r1 = self.stationary_source.r1

    def pixel_coordinates(self, vertices):
        # using the self.vertices, figure out the conversion to pixel
        # coordinates that facilitates deposition.
        # basically if ri are your stationary source grid then you need pixel
        # coordinates to obey (ri - r0) / r1
        return (vertices - self.r0) / (self.r1 - self.r0)

    def move_fluxes(self, vertices, fluxes, **kwargs):
        return deposit(self.fluxes, self.pixel_coordinates(vertices))


class UniformIlluminationMover(StationaryMover):
    """
    Basically get the area of each vertex
    """
    def __init__(self, flux=1, **kwargs):
        self.flux = flux  # counts per area

    def move_fluxes(self, vertices, fluxes, **kwargs):
        # flux in a box is based on total area of box
        x = vertices[:, :, 0]
        y = vertices[:, :, 1]
        # each box has (e.g. x) x[i,i] x[i,i+1], x[i+1,i], x[i+1,i+1]
        fluxes = 0.5 * np.abs((x[:-1, :-1] - x[1:, 1:]) *
                              (y[:-1, 1:] -  y[1:, :-1]) -
                              (x[:-1, 1:] -  x[1:, :-1]) *
                              (y[:-1, :-1] - y[1:, 1:]))
        fluxes *= self.flux
        return fluxes

class GaussianVerticesMover(StationaryMover):
    """
    Given a covariance matrix and means, move vertices by sampling gaussian.
    """

    def __init__(self, mu_x=0, mu_y=0, sigma_xx=1, sigma_yy=1, sigma_xy=0,
                 **kwargs):
        super(GaussianVerticesMover, self).__init__(**kwargs)
        self.mean = np.array([mu_x, mu_y])
        self.cov = np.array([[sigma_xx, sigma_xy],
                             [sigma_xy, sigma_yy]])
        # TODO: add check that cov is positive-semidefinite

    def move_vertices(self, vertices, fluxes, **kwargs):
        # realize perturbations from gaussian
        perts = np.random.multivariate_normal(self.mean, self.cov,
                size=vertices.shape[:2])
        return vertices + perts

class UniformGaussianMover(GaussianVerticesMover, UniformIlluminationMover):
    """
    multiple inheritence
    """
    def __init__(self, **kwargs):
        super(UniformGaussianMover, self).__init__(**kwargs)


class TreeringVerticesMover(StationaryMover):
    """
    Give a center, wavelength, amplitude, and phase.
    """

    def __init__(self, center=np.array([0,0]), wavelength=10,
                 amplitude=1, phase=0, **kwargs):
        super(TreeringVerticesMover, self).__init__(**kwargs)
        self.center = center
        self.wavelength = wavelength
        self.amplitude = amplitude
        self.phase = phase

    def move_vertices(self, vertices, fluxes, **kwargs):
        # convert vertices to radial coordinates about the center
        r = np.sqrt(np.sum(np.square(vertices - self.center), axis=2))
        theta = np.arctan2(vertices[:, :, 1] - self.center[1],
                           vertices[:, :, 0] - self.center[0])

        # displace along radius
        r += self.amplitude * np.sin(2 * np.pi * (r / self.wavelength +
                                                  self.phase))

        vertices = np.dstack((r * np.cos(theta),
                              r * np.sin(theta))) + self.center
        return vertices

class UniformTreeringMover(TreeringVerticesMover, UniformIlluminationMover):
    """
    multiple inheritence
    """
    def __init__(self, **kwargs):
        super(UniformTreeringMover, self).__init__(**kwargs)


if __name__ == '__main__':
    # show off treeringmover and plot it

    # make a gaussianmover, recover the stats doing something liek this:
    # np.cov((GaussianMover()(vertices) -
    # vertices).reshape((np.prod(vertices.shape[:2]), vertices.shape[2])).T)
    pass
