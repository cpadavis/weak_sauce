"""
Movers update flux AND vertices

TODO: Add method for adding movers
TODO: Make sure all move methods are now differences, instead of total moves
"""

import numpy as np

from weak_sauce.r2d import deposit

class Mover(object):
    """
    Class shell for movers. Currently super super simple. Just needs to have
    some kind of moving function that it calls.
    """

    def __init__(self, **kwargs):
        #print('M', kwargs)
        pass

    def move_vertices(self, vertices, fluxes, **kwargs):
        # returns dvertices
        raise NotImplementedError

    def move_fluxes(self, vertices, fluxes, **kwargs):
        # returns dfluxes
        raise NotImplementedError

    def move(self, source, **kwargs):
        # step order is to move vertices, then move fluxes
        source.vertices += self.move_vertices(source.vertices, source.fluxes,
                                              **kwargs)
        source.update_centroids()
        source.fluxes += self.move_fluxes(source.vertices, source.fluxes,
                                          **kwargs)

    def merge(self, other):
        print('Warning! You are automagically adding two movers together! This might be awesome, but it might also be catastrophic!')

        # so I didn't think this would work. If you just put in the methods,
        # you get a recursion error
        self_move_vertices = self.move_vertices
        other_move_vertices = other.move_vertices
        self_move_fluxes = self.move_fluxes
        other_move_fluxes = other.move_fluxes
        def merged_vertices(vertices, fluxes, **kwargs):
            dvertices_self = self_move_vertices(vertices, fluxes, **kwargs)
            dvertices_other = other_move_vertices(vertices, fluxes, **kwargs)
            return dvertices_self + dvertices_other
        def merged_fluxes(vertices, fluxes, **kwargs):
            dfluxes_self = self_move_fluxes(vertices, fluxes, **kwargs)
            dfluxes_other = other_move_fluxes(vertices, fluxes, **kwargs)
            return dfluxes_self + dfluxes_other

        self.move_vertices = merged_vertices
        self.move_fluxes = merged_fluxes

    def __call__(self, source, **kwargs):
        self.move(source, **kwargs)

    def __add__(self, other):
        self.merge(other)
        return self

class StationaryMover(Mover):
    """
    Don't move anything!
    """
    def __init__(self, **kwargs):
        #print('SM', kwargs)
        super(StationaryMover, self).__init__(**kwargs)

    def move_vertices(self, vertices, fluxes, **kwargs):
        return 0

    def move_fluxes(self, vertices, fluxes, **kwargs):
        return 0

class FixedIlluminationMover(StationaryMover):
    """
    Deposit some fixed but funny-shaped OTHER source
    """
    def __init__(self, stationary_source):
        super(FixedIlluminationMover, self).__init__(**kwargs)
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
        dfluxes = deposit(self.fluxes, self.pixel_coordinates(vertices))
        return dfluxes


class UniformIlluminationMover(StationaryMover):
    """
    Basically get the area of each vertex
    """
    def __init__(self, luminosity=1, **kwargs):
        #print('UIM', kwargs)
        super(UniformIlluminationMover, self).__init__(**kwargs)
        self.luminosity = luminosity  # counts per area

    def area(self, vertices):
        """
        unsigned area of quadrilateral given vertices
        A_i,j = ((x_i+1,j+1 - x_i,j) * (y_i,j+1 - y_i+1,j)
              -  (x_i,j+1 - x_i+1,j) * (y_i+1,j+1 - y_i,j)) / 2
        """
        x = vertices[:, :, 0]
        y = vertices[:, :, 1]

        A = ((x[:-1, :-1] - x[1:, 1:]) *
             (y[:-1, 1:] -  y[1:, :-1]) -
             (x[:-1, 1:] -  x[1:, :-1]) *
             (y[:-1, :-1] - y[1:, 1:])) * 0.5

        return A

    def move_fluxes(self, vertices, fluxes, **kwargs):
        # flux in a box is based on total area of box
        # each box has (e.g. x) x[i,i] x[i,i+1], x[i+1,i], x[i+1,i+1]
        dfluxes = self.luminosity * np.abs(self.area(vertices))
        return dfluxes

class SimpleVerticesMover(StationaryMover):
    """
    Move every vertex by some uniform amount
    """
    def __init__(self, mu_x=0, mu_y=0, **kwargs):
        # print('SVM', kwargs)
        super(SimpleVerticesMover, self).__init__(**kwargs)
        self.mean = np.array([mu_x, mu_y])

    def move_vertices(self, vertices, fluxes, **kwargs):
        dvertices = self.mean
        return dvertices


class GaussianVerticesMover(StationaryMover):
    """
    Given a covariance matrix and means, move vertices by sampling gaussian.
    """

    def __init__(self, mu_x=0, mu_y=0, sigma_xx=1, sigma_yy=1, sigma_xy=0,
                 **kwargs):
        # print('GVM', kwargs)
        super(GaussianVerticesMover, self).__init__(**kwargs)
        self.mean = np.array([mu_x, mu_y])
        self.cov = np.array([[sigma_xx, sigma_xy],
                             [sigma_xy, sigma_yy]])
        # TODO: add check that cov is positive-semidefinite

    def move_vertices(self, vertices, fluxes, **kwargs):
        # realize perturbations from gaussian
        dvertices = np.random.multivariate_normal(self.mean, self.cov,
                                                  size=vertices.shape[:2])
        return dvertices

class UniformGaussianMover(GaussianVerticesMover, UniformIlluminationMover):
    """
    multiple inheritence. This class perturbs the vertices and then deposits
    flux assuming uniform illumination
    """
    def __init__(self, mu_x=0, mu_y=0, sigma_xx=1, sigma_yy=1, sigma_xy=0,
                 luminosity=1, **kwargs):
        # print('UGM', kwargs)
        super(UniformGaussianMover, self).__init__(mu_x=mu_x, mu_y=mu_y,
              sigma_xx=sigma_xx, sigma_yy=sigma_yy, sigma_xy=sigma_xy,
              luminosity=luminosity, **kwargs)



class TreeringVerticesMover(StationaryMover):
    """
    Give a center, wavelength, amplitude, and phase.
    """

    def __init__(self, center=np.array([0,0]), wavelength=10,
                 amplitude=1, phase=0, **kwargs):
        # print('TVM', kwargs)
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
        dr = self.amplitude * np.sin(2 * np.pi * (r / self.wavelength +
                                                  self.phase))

        dvertices = np.dstack((dr * np.cos(theta),
                               dr * np.sin(theta)))
        return dvertices

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


"""
TODO: Put this in a unit test!

# The below code snippet tests that addition of mover classes works.
# You should see in SVM1 that the x components of the vertices are 0.1
# SVM2 increments y components by 0.1
# SVM12 should increment both x and y by 0.1

from weak_sauce.movers import SimpleVerticesMover as SVM
from weak_sauce.sources import Source
SVM0 = SVM()
SVM1 = SVM(mu_x = 0.1)
SVM2 = SVM(mu_y = 0.1)
SVM12 = SVM(mu_x = 0.1)
SVM12_2 = SVM(mu_x = 0.1)
SVM12.merge(SVM2)
SVM12_2 = SVM12_2 + SVM2
source0 = Source(num_x=2)
print(source0.vertices)
source = Source(num_x=2)
SVM0(source)
print(source.vertices - source0.vertices)
source = Source(num_x=2)
SVM1(source)
print(source.vertices - source0.vertices)
source = Source(num_x=2)
SVM2(source)
print(source.vertices - source0.vertices)
source = Source(num_x=2)
SVM12(source)
print(source.vertices - source0.vertices)
source = Source(num_x=2)
SVM12_2(source)
print(source.vertices - source0.vertices)

# test multiple inheritance
from weak_sauce.movers import UniformGaussianMover
UGM = UniformGaussianMover(luminosity=2, sigma_xx=0, sigma_yy=0)
print(UGM.luminosity)

from weak_sauce.sources import Source
source = Source(num_x=2)

print(source.vertices)
print(source.fluxes)
UGM(source)
print(source.vertices)
print(source.fluxes)
"""
