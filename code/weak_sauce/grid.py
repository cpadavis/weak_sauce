"""
In here goes the grid class.
"""

import numpy as np
from weak_sauce.sources import check_vertices

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

    def move_vertices(self, **kwargs):
        # some method that modifies self.vertices
        self.source.vertices += self.mover.move_vertices(
                self.source.vertices, self.source.fluxes, **kwargs)

    def move_fluxes(self, **kwargs):
        # some method that updates fluxes. This can depend on the vertices
        # and also on some external function
        # This basically calls self.source function and adds it to the
        # fluxes
        self.source.fluxes += self.mover.move_fluxes(
                self.source.vertices, self.source.fluxes, **kwargs)

    def step(self, **kwargs):
        # self.update_fluxes(**kwargs)
        # self.update_vertices(**kwargs)
        self.mover(self.source, **kwargs)

    def lnlike(self, **kwargs):
        raise NotImplementedError

    def fit(self, xtol=0.0001, ftol=0.0001, maxiter=10000, maxfun=None,verbose=False,
            **kwargs):
        """
        ala scipy.optimize:
        xtol : float, optional
            Relative error in xopt acceptable for convergence.
        ftol : number, optional
            Relative error in func(xopt) acceptable for convergence.
        maxiter : int, optional
            Maximum number of iterations to perform.
        maxfun : number, optional
            Maximum number of function evaluations to make.
            TODO: Currently not implimented
        """
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
        loss_history = []
        average_delta_param_history = []
        for it in xrange(maxiter):
            vertices_old = self.source.vertices.copy()
            fluxes_old = self.source.fluxes.copy()
            lnlike_old = np.abs(self.lnlike(**kwargs).copy())
            self.step(**kwargs)
            lnlike = np.abs(self.lnlike(**kwargs))
            if verbose: print (lnlike)
            loss_history.append(lnlike)

            delta_vx = np.sqrt(np.mean(np.square(self.source.vertices -
                vertices_old)[:,:,0]))
            delta_vy = np.sqrt(np.mean(np.square(self.source.vertices -
                vertices_old)[:,:,1]))
            delta_fluxes = np.sqrt(np.mean(np.square(self.source.fluxes -
                                                     fluxes_old)))
            deltas = np.array([delta_vx, delta_vy, delta_fluxes])
            if verbose: print (deltas)
            average_delta_param_history.append(deltas)
            # check changes
            if np.abs(lnlike - lnlike_old)/lnlike < ftol:
                print ('ftol reached')
                break
            if np.any(deltas < xtol):
                print ('xtol reached')
                break

    # wrap to the source object
    def evaluate_psf(self):
        # evaluate moments of fluxes image naievely.
        return self.source.evaluate_psf()

    def plot(self, **kwargs):
        return self.source.plot(**kwargs)

    def plot_real_grid(self, fig=None, ax=None):
        return self.source.plot_real_grid(fig=fig, ax=ax)

    def plot_pixel_grid(self, fig=None, ax=None):
        return self.source.plot_naieve_grid(fig=fig, ax=ax)

    def plot_vertices(self, fig=None, ax=None):
        return self.plot_vertices.plot_naieve_grid(fig=fig, ax=ax)
