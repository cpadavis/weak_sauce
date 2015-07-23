"""
In here goes the grid class.
"""

import numpy as np
from weak_sauce.sources import check_vertices
import pickle

class MoveableGrid(object):
    """A mapper between two sets of grids -- say pixel and focal plane.

    Vertices is an M x N x 2 array. As long as pixels don't cross (and they
    better not!), all boxes are defined as box[m, n] having vertices v[m, n],
    v[m + 1, n], v[m, n + 1], and v[m + 1, n + 1]. The 2 are for x and y.

    Terminology: real grids are the actual coordinates of the vertices.
                 ideal grids are from the indices of the fluxes
    """

    def __init__(self, *args, **kwargs):
        """
        Instantiates a MoveableGrid in one of two ways:

        1. __init__(pickle_name)
            pickle_name: a string representing the filename of some pickle
            created by MoveableGrid.saveto(pickle_name)
        2. __init__(source,mover)
            where source, mover have already been constructed.
        """

        if len(args) == 1:
            temp = pickle.load(open(args[0]))
            self.source = temp.source
            self.mover  = temp.mover
        elif len(args) == 2:
            #check ordering--mover should have move method
            assert hasattr(args[1],'move')
            self.source = args[0]
            self.mover  = args[1]
        else: raise Exception('__init__ takes 1 or 2 arguments!')

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

    def fit(self, xtol=1e-8, ftol=1e-6, maxiter=10000, step_size=None,
            maxfun=None, verbose=False, learning_rate_decay=0,
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

        learning_rate_decay : if step_size is specified, after every update, multiply step_size by (1 - learning_rate_decay)
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
        self.loss_history = [self.lnlike(**kwargs)]  # number, not object
        self.average_relative_delta_param_history = []
        for it in xrange(maxiter):
            if verbose: print (it)
            vertices_old = self.source.vertices.copy()
            fluxes_old = self.source.fluxes.copy()
            lnlike_old = self.loss_history[-1]
            self.step(step_size=step_size, **kwargs)
            lnlike = self.lnlike(**kwargs)
            if verbose: print (lnlike)
            self.loss_history.append(lnlike)

            delta_vx = np.sqrt(np.mean(np.square((self.source.vertices -
                vertices_old) / self.source.vertices)[:,:,0]))
            delta_vy = np.sqrt(np.mean(np.square((self.source.vertices -
                vertices_old) / self.source.vertices)[:,:,1]))
            delta_fluxes = np.sqrt(np.mean(np.square((self.source.fluxes -
                fluxes_old) / self.source.fluxes)))
            deltas = np.array([delta_vx, delta_vy, delta_fluxes])
            if verbose: print (deltas)
            self.average_relative_delta_param_history.append(deltas)


            if type(step_size) != type(None):
                step_size *= 1 - learning_rate_decay

            # check changes
            if np.abs((lnlike - lnlike_old) / lnlike) < ftol:
                print ('ftol reached')
                return
            if np.all(deltas < xtol):
                print ('xtol reached')
                return
        print('maxiter reached')
        return

    # wrap to the source object
    def evaluate_psf(self):
        # evaluate moments of fluxes image naievely.
        return self.source.evaluate_psf()

    def evaluate_sex(self,**kwargs):
        return self.source.evaluate_sex(**kwargs)

    def plot(self, **kwargs):
        return self.source.plot(**kwargs)

    def plot_real_grid(self, fig=None, ax=None):
        return self.source.plot_real_grid(fig=fig, ax=ax)

    def plot_pixel_grid(self, fig=None, ax=None):
        return self.source.plot_naieve_grid(fig=fig, ax=ax)

    def plot_vertices(self, fig=None, ax=None):
        return self.plot_vertices.plot_naieve_grid(fig=fig, ax=ax)

    def saveto(self,filename):
        file = open(filename,'w')
        pickle.dump(self,file)
        return
