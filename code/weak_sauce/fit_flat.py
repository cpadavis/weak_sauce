"""
Given a set of true fluxes, move vertices to try to match it.

TODO: Different types of movings
TODO: Convergence Criteria (to be put into a moveable grid
TODO: Keep track of relative changes?
TODO: FlatFitter should be a moveable grid. The routine to create a step for moving the vertices can be a Mover still.
"""


import numpy as np


from weak_sauce.movers import UniformIlluminationMover
from weak_sauce.grid import MoveableGrid

class FlatFitter(MoveableGrid):
    def __init__(self, source, true_fluxes, luminosity=1, step_size=1e-4,
                 **kwargs):
        mover = FlatMover(true_fluxes=true_fluxes, luminosity=luminosity,
                          step_size=step_size)
        super(FlatFitter, self).__init__(source=source, mover=mover, **kwargs)

    def lnlike(self):
        return self.mover.lnlike(self.source.vertices, self.source.fluxes)



class FlatMover(UniformIlluminationMover):
    """
    This Mover takes true_fluxes and uses them to tell the source how to move.
    move_fluxes method based on UniformIlluminationMover
    """
    def __init__(self, true_fluxes, luminosity=1, step_size=1e-4, **kwargs):
        super(FlatMover, self).__init__(luminosity=luminosity, **kwargs)
        self.true_fluxes = true_fluxes
        self.step_size = step_size

    def lnlike(self, vertices, fluxes):
        return -0.5 * np.sum(np.square(fluxes - self.true_fluxes))

    def derivatives(self, vertices, fluxes):
        """
        lnlike = -0.5 * sum_ij((A_ij - L * Ahat_ij)^2)
        """
        luminosity = np.abs(self.luminosity)
        x = vertices[:, :, 0]
        y = vertices[:, :, 1]

        # weight
        A = self.area(vertices)
        w_ij = (self.true_fluxes - luminosity * np.abs(A)) * np.sign(A) * luminosity

        # each displacement
        dA_ij__dx_ij = (y[:-1, 1:] - y[1:, :-1]) * 0.5
        dA_ij__dx_ip1jp1 = - (y[:-1, 1:] - y[1:, :-1]) * 0.5
        dA_ij__dx_ijp1 = - (y[:-1, :-1] - y[1:, 1:]) * 0.5
        dA_ij__dx_ip1j = (y[:-1, :-1] - y[1:, 1:]) * 0.5

        dA_ij__dy_ij = - (x[:-1, 1:] - x[1:, :-1]) * 0.5
        dA_ij__dy_ip1jp1 = (x[:-1, 1:] - x[1:, :-1]) * 0.5
        dA_ij__dy_ijp1 = (x[:-1, :-1] - x[1:, 1:]) * 0.5
        dA_ij__dy_ip1j = - (x[:-1, :-1] - x[1:, 1:]) * 0.5

        # now combine these so we can get in terms of dx_ij only
        """
        for i in xrange(Nx):
            for j in xrange(Ny):
                dLdx[i, j] += dA_ij__dx_ij[i, j]
                dLdx[i + 1, j + 1] += dA_ij__dx_ip1jp1[i, j]
                dLdx[i, j + 1] += dA_ij__dx_ijp1[i, j]
                dLdx[i + 1, j] += dA_ij__dx_ip1j[i, j]

        for alpha in xrange(Nx):
            for beta in xrange(Ny):
                dLdx[alpha, beta] += dA__ij__dx_ij[alpha, beta]
                dLdx[alpha, beta] += dA__ij__dx_ip1jp1[alpha - 1, beta - 1]
                dLdx[alpha, beta] += dA__ij__dx_ijp1[alpha, beta - 1]
                dLdx[alpha, beta] += dA__ij__dx_ip1j[alpha - 1, beta]
        """
        dLdx = np.zeros(x.shape)
        dLdx[:-1, :-1] += dA_ij__dx_ij * w_ij
        dLdx[1:, 1:] += dA_ij__dx_ip1jp1 * w_ij
        dLdx[:-1, 1:] += dA_ij__dx_ijp1 * w_ij
        dLdx[1:, :-1] += dA_ij__dx_ip1j * w_ij

        dLdy = np.zeros(y.shape)
        dLdy[:-1, :-1] += dA_ij__dy_ij * w_ij
        dLdy[1:, 1:] += dA_ij__dy_ip1jp1 * w_ij
        dLdy[:-1, 1:] += dA_ij__dy_ijp1 * w_ij
        dLdy[1:, :-1] += dA_ij__dy_ip1j * w_ij

        dLdvertices = np.dstack((dLdx, dLdy))

        # # TODO: This normalization by sum of true fluxes cannot be correct
        # dLdluminosity = np.sum((self.true_fluxes - luminosity * np.abs(A)) *
        #                        np.abs(A)) / np.sum(self.true_fluxes)

        return dLdvertices

    def move_vertices(self, vertices, fluxes, step_size=None,
                      **kwargs):

        if type(step_size) == type(None):
            step_size = self.step_size
        # get and apply derivatives
        dLdvertices = self.derivatives(vertices, fluxes)
        # apply derivatives
        dvertices = step_size * dLdvertices

        # TODO: This also doesn't work.
        #self.luminosity = self.luminosity + step_size * self.dLdluminosity

        # enforce that two edges (say the X edges) cannot move
        # enforce by translation and squish and compensation in luminosity
        # TODO: This didn't really work and actually was causing bugs. Shit.
        # # sort vertices to maintain ordering
        # v_flat = vertices.reshape(vertices.shape[0] * vertices.shape[1],
        #                           vertices.shape[2])
        # # non-trivial
        # vertices = np.sort(v_flat[np.argsort(v_flat[:,0])].reshape(
        #     vertices.shape), axis=1)

        return dvertices

    def move_fluxes(self, vertices, fluxes, **kwargs):
        # old move_fluxes adds a new flux with time (ie integrate over some time step)
        # but here we are not changing time, but the actual area of pixels
        # so we want the difference between the two
        updated_fluxes = super(FlatMover, self).move_fluxes(vertices, fluxes, **kwargs)
        return updated_fluxes - fluxes
