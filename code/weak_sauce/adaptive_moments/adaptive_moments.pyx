from __future__ import print_function, division
from math import atan2, cos, sin
import numpy as np
# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np
from libc.math cimport sqrt, ceil, floor, exp


pixel_to_arcsecond = 0.27

"""
Adaptation of GalSim's adaptive moments code, which is itself an adaptation of
some other code.

441 mus vs 338 mus for galsim python interface + cpp code (and I only cythonized ellipmom_1 !)
but: that is for a weaker epsilon of 1e-3. 1e-4 makes it take much longer (hits num_iter max)

TODO: why is it that I need that /2 in my ellipticity calculations to
recover roughly sextractor's ellipticity values? (for when I switch from
moments to ellipticities...) BUT not in e0!
"""

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

"A parameter for optimizing calculations of adaptive moments by\n"
"cutting off profiles. This parameter is used to decide how many\n"
"sigma^2 into the Gaussian adaptive moment to extend the moment\n"
"calculation, with the weight being defined as 0 beyond this point.\n"
"i.e., if max_moment_nsig2 is set to 25, then the Gaussian is\n"
"extended to (r^2/sigma^2)=25, with proper accounting for elliptical\n"
"geometry.  If this parameter is set to some very large number, then\n"
"the weight is never set to zero and the exponential function is\n"
"always called. Note: GalSim script devel/modules/test_mom_timing.py\n"
"was used to choose a value of 25 as being optimal, in that for the\n"
"cases that were tested, the speedups were typically factors of\n"
"several, but the results of moments and shear estimation were\n"
"changed by <10^-5.  Not all possible cases were checked, and so for\n"
"use of this code for unusual cases, we recommend that users check\n"
"that this value does not affect accuracy, and/or set it to some\n"
"large value to completely disable this optimization.\n"
cdef double MAX_MOMENT_NSIG2 = 25.

# cpdef NOT cdef !!
cpdef np.ndarray[DTYPE_t, ndim=1] find_ellipmom_1(
        np.ndarray[DTYPE_t, ndim=2] data,
        double Mx, double My,
        double Mxx, double Mxy, double Myy):
    """C routine to obtain estimates of corrections to Mxx etc. Taken from
    galsim."""


    cdef DTYPE_t A = 0
    cdef DTYPE_t Bx = 0
    cdef DTYPE_t By = 0
    cdef DTYPE_t Cxx = 0
    cdef DTYPE_t Cxy = 0
    cdef DTYPE_t Cyy = 0
    cdef DTYPE_t Cxxx = 0
    cdef DTYPE_t Cxxy = 0
    cdef DTYPE_t Cxyy = 0
    cdef DTYPE_t Cyyy = 0
    cdef DTYPE_t rho4w = 0

    cdef np.ndarray[DTYPE_t, ndim=1] return_array

    cdef int xmin = 0
    cdef int ymin = 0
    cdef int ymax = data.shape[0]
    cdef int xmax = data.shape[1]

    cdef double detM = Mxx * Myy - Mxy * Mxy
    if (detM <= 0) + (Mxx <= 0) + (Myy <= 0):
        print("Error: non positive definite adaptive moments!\n")

    cdef double Minv_xx = Myy / detM
    cdef double TwoMinv_xy = -Mxy / detM * 2.0
    cdef double Minv_yy = Mxx / detM
    cdef double Inv2Minv_xx = 0.5 / Minv_xx  # Will be useful later...

    # rho2 = Minv_xx(x-Mx)^2 + 2Minv_xy(x-Mx)(y-My) + Minv_yy(y-My)^2
    # The minimum/maximum y that have a solution rho2 = max_moment_nsig2 is at:
    #   2*Minv_xx*(x-Mx) + 2Minv_xy(y-My) = 0
    # rho2 = Minv_xx (Minv_xy(y-My)/Minv_xx)^2
    #           - 2Minv_xy(Minv_xy(y-My)/Minv_xx)(y-My)
    #           + Minv_yy(y-My)^2
    #      = (Minv_xy^2/Minv_xx - 2Minv_xy^2/Minv_xx + Minv_yy) (y-My)^2
    #      = (Minv_xx Minv_yy - Minv_xy^2)/Minv_xx (y-My)^2
    #      = (1/detM) / Minv_xx (y-My)^2
    #      = (1/Myy) (y-My)^2
    #
    # we are finding the limits for the iy values and then the ix values.
    cdef double y_My = sqrt(MAX_MOMENT_NSIG2 * Myy)
    # nan check!
    if y_My != y_My:
        return_array = \
            np.array([A, Bx, By, Cxx, Cxy, Cyy, rho4w,
                      Cxxx, Cxxy, Cxyy, Cyyy
                      ], dtype=DTYPE)
        return return_array

    cdef double y1 = -y_My + My
    cdef double y2 = y_My + My

    # stay within image bounds
    cdef int iy1 = max(<int>(ceil(y1)), ymin)
    cdef int iy2 = min(<int>(floor(y2)), ymax)
    cdef int y

    if iy1 > iy2:
        print('iy1 > iy2', y1, ymin, y2, ymax, iy1, iy2)

    cdef double a, b, c, d, sqrtd, inv2a, x1, x2, x_Mx, \
        Minv_xx__x_Mx__x_Mx, rho2, intensity, TwoMinv_xy__y_My, Minv_yy__y_My__y_My
    cdef int ix1, ix2, x

    for y in xrange(iy1, iy2):

        y_My = float(y) - My
        TwoMinv_xy__y_My = TwoMinv_xy * y_My
        Minv_yy__y_My__y_My = Minv_yy * y_My ** 2

        # Now for a particular value of y, we want to find the min/max x that satisfy
        # rho2 < max_moment_nsig2.
        #
        # 0 = Minv_xx(x-Mx)^2 + 2Minv_xy(x-Mx)(y-My) + Minv_yy(y-My)^2 - max_moment_nsig2
        # Simple quadratic formula:

        a = Minv_xx
        b = TwoMinv_xy__y_My
        c = Minv_yy__y_My__y_My - MAX_MOMENT_NSIG2
        d = b * b - 4 * a * c
        sqrtd = sqrt(d)
        inv2a = Inv2Minv_xx
        x1 = inv2a * (-b - sqrtd) + Mx
        x2 = inv2a * (-b + sqrtd) + Mx

        # stay within image bounds
        ix1 = max(<int>(ceil(x1)), xmin)
        ix2 = min(<int>(floor(x2)), xmax)
        # in the following two cases, ask if we somehow wanted to find
        # pixels outside the image
        if (ix1 > xmax) * (ix2 == xmax):
            continue
        elif (ix1 == xmin) * (ix2 < xmin):
            continue
        elif ix1 > ix2:
            # print('ix1 > ix2', y, x1, xmin, x2, xmax, ix1, ix2)
            # usually what happens is you want to take only one pixel and you
            # end up due to the ceil and floor funcs with e.g. 15, 14 instead
            # of 14, 15
            # ix1, ix2 = ix2, ix1
            # ix1 = max(ix1, xmin)
            # ix2 = min(ix2, xmax)
            continue

        for x in xrange(ix1, ix2):

            x_Mx = float(x) - Mx

            # Compute displacement from weight centroid, then get elliptical
            # radius and weight.
            Minv_xx__x_Mx__x_Mx = Minv_xx * x_Mx ** 2
            rho2 = Minv_yy__y_My__y_My + \
                TwoMinv_xy__y_My * x_Mx + \
                Minv_xx__x_Mx__x_Mx

            # this shouldn't happen by construction
            if (rho2 > MAX_MOMENT_NSIG2 + 1e8):
                print('rho2 > max_moment_nsig2 !')
                continue

            intensity = exp(-0.5 * rho2) * data[y, x]  # y,x order!

            A += intensity
            Bx += intensity * x_Mx
            By += intensity * y_My
            Cxx += intensity * x_Mx ** 2
            Cxy += intensity * x_Mx * y_My
            Cyy += intensity * y_My ** 2
            Cxxx += intensity * x_Mx ** 3
            Cxxy += intensity * x_Mx ** 2 * y_My
            Cxyy += intensity * x_Mx * y_My ** 2
            Cyyy += intensity * y_My ** 3
            rho4w += intensity * rho2 * rho2

    return_array = \
        np.array([A, Bx, By, Cxx, Cxy, Cyy, rho4w,
                  Cxxx, Cxxy, Cxyy, Cyyy
                  ], dtype=DTYPE)
    return return_array

def adaptive_moments(data,
    epsilon = 1e-6,
    convergence_factor = 1.0,
    guess_sig = 3.0,
    bound_correct_wt = 0.25,  # Maximum shift in centroids and sigma between
                             # iterations for adaptive moments.
    guess_centroid = None
    num_iter = 0,
    num_iter_max = 100,
    **kwargs
    ):

    # Set Amp = -1000 as initial value just in case the while() block below is
    # never triggered; in this case we have at least *something* defined to
    # divide by, and for which the output will fairly clearly be junk.
    Amp = -1000.
    if type(guess_centroid) == None:
        Mx = data.shape[1] / 2.
        My = data.shape[0] / 2.
    else:
        Mx = guess_centroid[0]
        My = guess_centroid[1]
    Mxx = guess_sig ** 2
    Mxy = 0
    Myy = guess_sig ** 2

    # Iterate until we converge
    while (convergence_factor > epsilon) * (num_iter < num_iter_max):
        # print(Mx, My, Mxx, Mxy, Myy, num_iter)

        # Get moments
        Amp, Bx, By, Cxx, Cxy, Cyy, rho4, \
            Cxxx, Cxxy, Cxyy, Cyyy \
            = find_ellipmom_1(data, Mx, My, Mxx, Mxy, Myy)
        # print(num_iter)
        # print(Amp, Bx / Amp, By / Amp, Cxx / Amp, Cxy / Amp, Cyy / Amp, rho4 / Amp)
        # print(Mx, My, Mxx, Mxy, Myy)
        # Compute configuration of the weight function
        two_psi = atan2(2 * Mxy, Mxx - Myy)
        semi_a2 = 0.5 * ((Mxx + Myy) + (Mxx - Myy) * cos(two_psi)) + \
                         Mxy * sin(two_psi)
        semi_b2 = Mxx + Myy - semi_a2

        if semi_b2 <= 0:
            print("Error: non positive-definite weight in find_ellipmom_2.\n")

        shiftscale = sqrt(semi_b2)
        if num_iter == 0:
            shiftscale0 = shiftscale

        # Now compute changes to Mx, etc
        dx = 2. * Bx / (Amp * shiftscale)
        dy = 2. * By / (Amp * shiftscale)
        dxx = 4. * (Cxx / Amp - 0.5 * Mxx) / semi_b2
        dxy = 4. * (Cxy / Amp - 0.5 * Mxy) / semi_b2
        dyy = 4. * (Cyy / Amp - 0.5 * Myy) / semi_b2

        if (dx     >  bound_correct_wt): dx     =  bound_correct_wt
        if (dx     < -bound_correct_wt): dx     = -bound_correct_wt
        if (dy     >  bound_correct_wt): dy     =  bound_correct_wt
        if (dy     < -bound_correct_wt): dy     = -bound_correct_wt
        if (dxx    >  bound_correct_wt): dxx    =  bound_correct_wt
        if (dxx    < -bound_correct_wt): dxx    = -bound_correct_wt
        if (dxy    >  bound_correct_wt): dxy    =  bound_correct_wt
        if (dxy    < -bound_correct_wt): dxy    = -bound_correct_wt
        if (dyy    >  bound_correct_wt): dyy    =  bound_correct_wt
        if (dyy    < -bound_correct_wt): dyy    = -bound_correct_wt

        # Convergence tests
        convergence_factor = abs(dx) ** 2
        if (abs(dy) > convergence_factor):
            convergence_factor = abs(dy) ** 2
        if (abs(dxx) > convergence_factor):
            convergence_factor = abs(dxx)
        if (abs(dxy) > convergence_factor):
            convergence_factor = abs(dxy)
        if (abs(dyy) > convergence_factor):
            convergence_factor = abs(dyy)
        convergence_factor = np.sqrt(convergence_factor)
        if (shiftscale < shiftscale0):
            convergence_factor *= shiftscale0 / shiftscale

        # Now update moments
        Mx += dx * shiftscale
        My += dy * shiftscale
        Mxx += dxx * semi_b2
        Mxy += dxy * semi_b2
        Myy += dyy * semi_b2

        num_iter += 1

    # we made it! do a final calculation
    Amp, Bx, By, Cxx, Cxy, Cyy, rho4, \
        Cxxx, Cxxy, Cxyy, Cyyy \
        = find_ellipmom_1(data, Mx, My, Mxx, Mxy, Myy)

    A = Amp
    rho4 /= Amp
    x2 = Cxx / Amp
    xy = Cxy / Amp
    y2 = Cyy / Amp
    x3 = Cxxx / Amp
    x2y = Cxxy / Amp
    xy2 = Cxyy / Amp
    y3 = Cyyy / Amp

    return Mx, My, Mxx, Mxy, Myy, A, rho4, x2, xy, y2, x3, x2y, xy2, y3


cpdef double centered_moment(
        np.ndarray[DTYPE_t, ndim=2] data,
        int p, int q,
        double Mx, double My,
        double Mxx, double Mxy, double Myy):

    """Calculate centered moments

    Parameters
    ----------
    data : array
        2d image array

    p, q : int
        The moments we are interested in. p,q for x,y

    Mx, My : floats
        The centroids of the image.

    Mxx, Mxy, Myy : floats
        The second moment weight matrix. See the notes below on the extra
        factor of two compared to the actual centered moment.

    Returns
    -------
    mu_pq : float
        The centered adaptive moment.

    Notes
    -----
    Because of the 2 in eq (5) of Hirata et al 2004, the second moment Ms in
    the input are /twice/ the value of the centered moments I would obtain
    otherwise. When you measure ellipticities it generally doesn't matter since
    you normalize by Mxx + Myy, thus cancelling the factor of two. In my
    analysis, however, I have kept these separate so this is good to
    understand.

    ie Mxx = 2 <x^2>, Mxy = 2 <xy>, Myy = 2 <y^2>
    but Mx = <x>, My = <y>

    consequently the T/2 in (4) = e0 = <x^2> + <y^2>

    """

    cdef double A = 0
    cdef double mu_pq = 0

    cdef int xmin = 0
    cdef int ymin = 0
    cdef int ymax = data.shape[0]
    cdef int xmax = data.shape[1]

    cdef double detM = Mxx * Myy - Mxy * Mxy
    if (detM <= 0) + (Mxx <= 0) + (Myy <= 0):
        print("Error: non positive definite adaptive moments!\n")

    cdef double Minv_xx = Myy / detM
    cdef double TwoMinv_xy = -Mxy / detM * 2.0
    cdef double Minv_yy = Mxx / detM
    cdef double Inv2Minv_xx = 0.5 / Minv_xx  # Will be useful later...

    # rho2 = Minv_xx(x-Mx)^2 + 2Minv_xy(x-Mx)(y-My) + Minv_yy(y-My)^2
    # The minimum/maximum y that have a solution rho2 = max_moment_nsig2 is at:
    #   2*Minv_xx*(x-Mx) + 2Minv_xy(y-My) = 0
    # rho2 = Minv_xx (Minv_xy(y-My)/Minv_xx)^2
    #           - 2Minv_xy(Minv_xy(y-My)/Minv_xx)(y-My)
    #           + Minv_yy(y-My)^2
    #      = (Minv_xy^2/Minv_xx - 2Minv_xy^2/Minv_xx + Minv_yy) (y-My)^2
    #      = (Minv_xx Minv_yy - Minv_xy^2)/Minv_xx (y-My)^2
    #      = (1/detM) / Minv_xx (y-My)^2
    #      = (1/Myy) (y-My)^2
    #
    # we are finding the limits for the iy values and then the ix values.
    cdef double y_My = sqrt(MAX_MOMENT_NSIG2 * Myy)
    # nan check!
    if y_My != y_My:
        return 0

    cdef double y1 = -y_My + My
    cdef double y2 = y_My + My

    # stay within image bounds
    cdef int iy1 = max(<int>(ceil(y1)), ymin)
    cdef int iy2 = min(<int>(floor(y2)), ymax)
    cdef int y
    if iy1 > iy2:
        print('iy1 > iy2', y1, ymin, y2, ymax, iy1, iy2)

    cdef double a, b, c, d, sqrtd, inv2a, x1, x2, x_Mx, \
        Minv_xx__x_Mx__x_Mx, rho2, intensity, \
        TwoMinv_xy__y_My, Minv_yy__y_My__y_My, \
        y_My_q
    cdef int ix1, ix2, x

    for y in xrange(iy1, iy2):

        y_My = float(y) - My
        y_My_q = y_My ** q

        TwoMinv_xy__y_My = TwoMinv_xy * y_My
        Minv_yy__y_My__y_My = Minv_yy * y_My ** 2

        # Now for a particular value of y, we want to find the min/max x that satisfy
        # rho2 < max_moment_nsig2.
        #
        # 0 = Minv_xx(x-Mx)^2 + 2Minv_xy(x-Mx)(y-My) + Minv_yy(y-My)^2 - max_moment_nsig2
        # Simple quadratic formula:

        a = Minv_xx
        b = TwoMinv_xy__y_My
        c = Minv_yy__y_My__y_My - MAX_MOMENT_NSIG2
        d = b * b - 4 * a * c
        sqrtd = sqrt(d)
        inv2a = Inv2Minv_xx
        x1 = inv2a * (-b - sqrtd) + Mx
        x2 = inv2a * (-b + sqrtd) + Mx

        # stay within image bounds
        ix1 = max(<int>(ceil(x1)), xmin)
        ix2 = min(<int>(floor(x2)), xmax)
        # in the following two cases, ask if we somehow wanted to find
        # pixels outside the image
        if (ix1 > xmax) * (ix2 == xmax):
            continue
        elif (ix1 == xmin) * (ix2 < xmin):
            continue
        elif ix1 > ix2:
            # print('ix1 > ix2', y, x1, xmin, x2, xmax, ix1, ix2)
            # usually what happens is you want to take only one pixel and you
            # end up due to the ceil and floor funcs with e.g. 15, 14 instead
            # of 14, 15
            # ix1, ix2 = ix2, ix1
            # ix1 = max(ix1, xmin)
            # ix2 = min(ix2, xmax)
            continue

        for x in xrange(ix1, ix2):

            x_Mx = float(x) - Mx

            # Compute displacement from weight centroid, then get elliptical
            # radius and weight.
            Minv_xx__x_Mx__x_Mx = Minv_xx * x_Mx ** 2
            rho2 = Minv_yy__y_My__y_My + \
                TwoMinv_xy__y_My * x_Mx + \
                Minv_xx__x_Mx__x_Mx

            # this shouldn't happen by construction
            if (rho2 > MAX_MOMENT_NSIG2 + 1e8):
                print('rho2 > max_moment_nsig2 !')
                continue

            intensity = exp(-0.5 * rho2) * data[y, x]  # y,x order!

            A += intensity

            if (p == 0) + (q == 0):
                if p == 0:
                    mu_pq += intensity * y_My_q
                else:
                    mu_pq += intensity * x_Mx ** p
            else:
                mu_pq += intensity * x_Mx ** p * y_My_q

    return mu_pq / A

def second_moment_to_ellipticity(x2, y2, xy, **args):

    """Take moments and convert to unnormalized ellipticity basis

    Parameters
    ----------
    x2, y2, xy : array
        Array of second moments, e.g. fits_data['X2WIN_IMAGE'] .

    Returns
    -------
    e0, e0prime, e1, e2 : array
        Arrays converted to unnormalized ellipticity basis.
        e0prime is an alternative definition for the spin 0 second moment.

    Notes
    -----
    If you want the normalized (and unitless) ellipticity:
        e1 -> e1 / e0
        e2 -> e2 / e0

    Units are arcsec ** 2 provided the moments are in pixel ** 2.

    Roughly, FWHM ~ sqrt(8 ln 2 e0).

    References
    ----------
    see http://des-docdb.fnal.gov:8080/cgi-bin/RetrieveFile?docid=353
    though my convention is different.

    """

    e0 = (x2 + y2) * pixel_to_arcsecond ** 2
    e0prime = (x2 + y2 + 2 * np.sqrt(x2 * y2 - xy ** 2)) * pixel_to_arcsecond ** 2
    e1 = (x2 - y2) * pixel_to_arcsecond ** 2
    e2 = (2 * xy) * pixel_to_arcsecond ** 2

    return e0, e0prime, e1, e2


def third_moments_to_octupoles(x3, x2y, xy2, y3, **args):

    """Take moments and convert to unnormalized octupole basis

    Parameters
    ----------
    x3, x2y, xy2, y3 : array
        Array of third moments.

    Returns
    -------
    zeta1, zeta2, delta1, delta2 : array
        Arrays converted to unnormalized octupole basis.
        zeta is spin-1, delta is spin-3 (See Okura 2008)
        these are F and G (roughly and modulo normalization factors)

    Notes
    -----

    Units are arcsec ** 3 provided the moments are in pixel ** 3.

    """

    zeta1 = (x3 + xy2) * pixel_to_arcsecond ** 3
    zeta2 = (y3 + x2y) * pixel_to_arcsecond ** 3
    delta1 = (x3 - 3 * xy2) * pixel_to_arcsecond ** 3
    delta2 = (-y3 + 3 * x2y) * pixel_to_arcsecond ** 3

    return zeta1, zeta2, delta1, delta2


def ellipticity_to_whisker(e1, e2, spin=2., power=2., **args):

    """Take unnormalized ellipticities and convert to relevant whisker
    parameters

    Parameters
    ----------
    e1, e2 : array
        Array of unnormalized ellipticity vectors.

    Returns
    -------
    u, v : array
        Whisker parameters in cartesian basis (for plotting with matplotlib
        quiver).

    w, phi : array
        Whisker length and position angle. (See reference below.)

    Notes
    -----
    If your unnormalized ellipticities are in quantity x$^{2}$ (e.g.
    pixels$^{2}$), then all these return quantities are in x (except for phi,
    which is in radians).

    References
    ----------
    see http://des-docdb.fnal.gov:8080/cgi-bin/RetrieveFile?docid=353

    """

    # get magnitude and direction
    w = np.sqrt(np.square(e1) + np.square(e2)) ** (1 / power)
    phi = np.arctan2(e2, e1) / spin
    # convert to cartesian
    u = w * np.cos(phi)
    v = w * np.sin(phi)

    return u, v, w, phi


def second_moment_variance_to_ellipticity_variance(var_x2, var_y2, var_xy,
                                                   **args):

    """Convert variance in moments to ellipticity variances

    Parameters
    ----------
    var_x2, var_y2, var_xy : array
        Array of second moments variances

    Returns
    -------
    var_e0, var_e1, var_e2 : array
        Arrays converted to unnormalized ellipticity basis.

    """

    alpha = pixel_to_arcsecond ** 2  # pixel to arcsec
    var_e0 = alpha ** 2 * (var_x2 + var_y2)
    var_e1 = alpha ** 2 * (var_x2 + var_y2)
    var_e2 = alpha ** 2 * 4 * var_xy

    return var_e0, var_e1, var_e2

def third_moment_variance_to_octupole_variance(
        var_x3, var_x2y, var_xy2, var_y3, **args):

    """Convert variance in moments to ellipticity variances assuming no
    covariance.

    Parameters
    ----------
    var_x3, etc : array
        Array of third moments variances

    Returns
    -------
    var_delta1 etc : array
        Arrays converted to unnormalized octupole basis.

    """

    alpha = pixel_to_arcsecond ** 3  # pixel to arcsec
    var_zeta1 = alpha ** 2 * (var_x3 + var_xy2)
    var_zeta2 = alpha ** 2 * (var_y3 + var_x2y)

    var_delta1 = alpha ** 2 * (var_x3 + 9 * var_xy2)
    var_delta2 = alpha ** 2 * (var_y3 + 9 * var_x2y)

    return var_zeta1, var_zeta2, var_delta1, var_delta2


def ellipticity_variance_to_whisker_variance(e1, e2, var_e1, var_e2, **args):

    """Convert error in ellipticity to error in cartesian whisker parameters
    (for later wedge plotting).

    Parameters
    ----------
    e1, e2 : array
        Array of unnormalized ellipticity vectors.

    var_e1, var_e2 : array or float
        Variance in measurement of e1 or e2, either one value for all, or one
        value per object.

    Returns
    -------
    var_u, var_v : array
        Variance in the u and v components of the whisker.

    """

    # get relevant whisker parameters
    u, v, w, phi = ellipticity_to_whisker(e1, e2)

    # thanks to the power of WolframAlpha:
    dude1 = 0.5 * w ** -4 * (e2 * v + e1 * u)
    dude2 = 0.5 * w ** -4 * (e2 * u - e1 * v)
    dvde1 = 0.5 * w ** -4 * (-e2 * u + e1 * v)
    dvde2 = 0.5 * w ** -4 * (e2 * v + e1 * u)

    # errors are added in quadrature
    var_u = dude1 ** 2 * var_e1 + dude2 ** 2 * var_e2
    var_v = dvde1 ** 2 * var_e1 + dvde2 ** 2 * var_e2

    return var_u, var_v

def convert_moments(data, **args):

    """Looks through data and converts all relevant parameters and updates into
    the data object.

    Parameters
    ----------
    data : recarray or dictionary
        Set of data with the moments

    Returns
    -------
    poles : dictionary
        Dictionary of results

    """

    poles = data.copy()
    if ('x2' in data) * ('y2' in data) * ('xy' in data):
        e0, e0prime, e1, e2 = second_moment_to_ellipticity(**data)
        poles.update(dict(e0=e0, e0prime=e0prime, e1=e1, e2=e2))

    if ('e1' in poles) * ('e2' in poles):
        w1, w2, w, phi = ellipticity_to_whisker(**poles)
        poles.update(dict(w1=w1, w2=w2, w=w, phi=phi))

    if ('x3' in data) * ('x2y' in data) * ('xy2' in data) * ('y3' in data):
        zeta1, zeta2, delta1, delta2 = third_moments_to_octupoles(**data)
        wd1, wd2 = ellipticity_to_whisker(delta1, delta2, spin=3, power=3)[:2]
        poles.update(dict(zeta1=zeta1, zeta2=zeta2,
                          delta1=delta1, delta2=delta2))

    if ('delta1' in poles) * ('delta2' in poles):
        wd1, wd2 = ellipticity_to_whisker(poles['delta1'], poles['delta2'],
                                          spin=3, power=3)[:2]
        poles.update(dict(wd1=wd1, wd2=wd2))


    if ('x4' in data) * ('x2y2' in data) * ('y4' in data):
        xi = data['x4'] + 2 * data['x2y2'] + data['y4']
        poles.update(dict(xi=xi))


    # Variances
    if (('var_x2' in data) * ('var_y2' in data) * ('var_xy' in data) *
            ('var_e0' not in data) * ('var_e1' not in data) *
            ('var_e2' not in data)):
        var_e0, var_e1, var_e2 = \
            second_moment_variance_to_ellipticity_variance(**data)
        poles.update(dict(var_e0=var_e0, var_e1=var_e1, var_e2=var_e2))

    if (('e1' in poles) * ('e2' in poles) *
            ('var_e1' in poles) * ('var_e2' in poles) *
            ('var_w1' not in poles) * ('var_w2' not in poles)):
        var_w1, var_w2 = ellipticity_variance_to_whisker_variance(**poles)
        poles.update(dict(var_w1=var_w1, var_w2=var_w2))

    if (('var_x3' in data) * ('var_x2y' in data) *
            ('var_xy2' in data) * ('var_y3' in data) *
            ('var_zeta1' not in data) * ('var_zeta2' not in data) *
            ('var_delta1' not in data) * ('var_delta2' not in data)):
        var_zeta1, var_zeta2, var_delta1, var_delta2 = \
            third_moment_variance_to_octupole_variance(**data)
        poles.update(dict(var_zeta1=var_zeta1, var_zeta2=var_zeta2,
                          var_delta1=var_delta1, var_delta2=var_delta2,))

    return poles
