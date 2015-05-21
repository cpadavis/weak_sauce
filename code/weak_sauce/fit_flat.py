"""
Given a set of true fluxes, move vertices to try to match it.

TODO: Different types of movings
TODO: Convergence Criteria (to be put into a moveable grid
TODO: Keep track of relative changes?
"""


import numpy as np


from weak_sauce.movers import Mover
from weak_sauce.grid import MoveableGrid

class FlatFitter(Mover):
    def __init__(self, true_fluxes, luminosity=1,tychoX=.1,tychoY=.1):
        self.true_fluxes = true_fluxes
        self.luminosity = luminosity
        self.dLdluminosity = 0
        self.firstStepTaken = False
        self.tychonoffSigmaX = tychoX
        self.tychonoffSigmaY = tychoY
        self.verbose = False

    def lnlike(self, vertices, fluxes):
        if self.tychonoffSigmaX != 0 and self.firstStepTaken :
            return -0.5 * np.sum(np.square(fluxes - self.true_fluxes)) \
            -0.5 * np.sum(np.square(vertices[:,:,0] - self.original_vertices[:,:,0]))/(self.tychonoffSigmaX**2) \
            -0.5 * np.sum(np.square(vertices[:,:,1] - self.original_vertices[:,:,1]))/(self.tychonoffSigmaY**2)
        else:
            return -0.5 * np.sum(np.square(fluxes - self.true_fluxes))

    def derivatives(self, vertices, fluxes):
        """
        lnlike = -0.5 * sum_ij((A_ij - L * Ahat_ij)^2)
        """
        luminosity = np.abs(self.luminosity)
        x = vertices[:, :, 0]
        y = vertices[:, :, 1]

        # weight
        A = self.area(x, y)
        w_ij = (self.true_fluxes - luminosity * np.abs(A)) * np.sign(A) * luminosity

        # each displacement
        dA_ij__dx_ij = (y[:-1, 1:] -  y[1:, :-1]) * 0.5
        dA_ij__dx_ip1jp1 = - (y[:-1, 1:] -  y[1:, :-1]) * 0.5
        dA_ij__dx_ijp1 = - (y[:-1, :-1] - y[1:, 1:]) * 0.5
        dA_ij__dx_ip1j = (y[:-1, :-1] - y[1:, 1:]) * 0.5

        dA_ij__dy_ij = - (x[:-1, 1:] -  x[1:, :-1]) * 0.5
        dA_ij__dy_ip1jp1 = (x[:-1, 1:] -  x[1:, :-1]) * 0.5
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
        tychonoff_adjustmentX = -0.5*np.square(self.original_vertices[:,:,0]-x)/(self.tychonoffSigmaX**2)
        if self.verbose:
            print 'x adjustment is: ', tychonoff_adjustmentX
        dLdx += tychonoff_adjustmentX

        dLdy = np.zeros(y.shape)
        dLdy[:-1, :-1] += dA_ij__dy_ij * w_ij
        dLdy[1:, 1:] += dA_ij__dy_ip1jp1 * w_ij
        dLdy[:-1, 1:] += dA_ij__dy_ijp1 * w_ij
        dLdy[1:, :-1] += dA_ij__dy_ip1j * w_ij
        tychonoff_adjustmentY = -0.5*np.square(self.original_vertices[:,:,1]-y)/(self.tychonoffSigmaY**2)
        dLdy += tychonoff_adjustmentY

        dLdvertices = np.dstack((dLdx, dLdy))

        # TODO: This normalization by sum of true fluxes cannot be correct
        self.dLdluminosity = np.sum((self.true_fluxes - luminosity * np.abs(A)) *
                               np.abs(A)) / np.sum(self.true_fluxes)

        return dLdvertices


    def area(self, x, y):
        """
        A_i,j = ((x_i+1,j+1 - x_i,j) * (y_i,j+1 - y_i+1,j)
              -  (x_i,j+1 - x_i+1,j) * (y_i+1,j+1 - y_i,j)) / 2
        """
        A = ((x[:-1, :-1] - x[1:, 1:]) *
             (y[:-1, 1:] -  y[1:, :-1]) -
             (x[:-1, 1:] -  x[1:, :-1]) *
             (y[:-1, :-1] - y[1:, 1:])) * 0.5
        return A

    def move_fluxes(self, vertices, fluxes, **kwargs):
        # flux in a box is based on total area of box
        x = vertices[:, :, 0]
        y = vertices[:, :, 1]
        # each box has (e.g. x) x[i,i] x[i,i+1], x[i+1,i], x[i+1,i+1]
        fluxes = self.luminosity * np.abs(self.area(x, y))
        return fluxes

    def move_vertices(self, vertices, fluxes, step_size,
                      **kwargs):
        # TODO: incorporate several different parameter update modes
        """
        if update == 'sgd':
          dx = -learning_rate * grads[p]
        elif update == 'momentum':
          if not p in self.step_cache:
            self.step_cache[p] = np.zeros(grads[p].shape)
          dx = np.zeros_like(grads[p]) # you can remove this after
          dx = momentum * self.step_cache[p] - learning_rate * grads[p]
          self.step_cache[p] = dx
        elif update == 'rmsprop':
          decay_rate = 0.99 # you could also make this an option
          if not p in self.step_cache:
            self.step_cache[p] = np.zeros(grads[p].shape)
          dx = np.zeros_like(grads[p]) # you can remove this after
          self.step_cache[p] = self.step_cache[p] * decay_rate + (1.0 - decay_rate) * grads[p] ** 2
          dx = -(learning_rate * grads[p]) / np.sqrt(self.step_cache[p] + 1e-8)
        """
        # get derivatives
        if not self.firstStepTaken:
            if self.verbose:
                print 'this is the first step'
            self.original_vertices = vertices.copy()
            self.firstStepTaken = True
        dLdvertices = self.derivatives(vertices, fluxes)
        # apply derivatives
        vertices = vertices + step_size * dLdvertices
        if self.verbose:
            print vertices
            print self.original_vertices
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

        return vertices
