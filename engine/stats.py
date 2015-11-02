import numpy as np  


def median_absolute_deviation(a, axis=None):
    """Compute the median absolute deviation

    Returns the median absolute deviation  of the array elements.  The MAD is
    defined as median(|a-median(a)|).

    Parameters
    ----------
    a : array_like
        Input array or object that can be converted to an array.
    axis : int, optional
        Axis along which the medians are computed. The default (axis=None)
        is to compute the median along a flattened version of the array.

    Returns
    -------
    median_absolute_deviation : ndarray
        A new array holding the result. If the input contains
        integers, or floats of smaller precision than 64, then the output
        data-type is float64.  Otherwise, the output data-type is the same
        as that of the input.

    Examples
    --------

    This will generate random variates from a Gaussian distribution and return
    the median absolute deviation for that distribution::

        >>> from astropy.stats import median_aboslute_deviation
        >>> from numpy.random import randn
        >>> randvar = randn(10000)
        >>> mad = median_absolute_deviation(randvar)

    See Also
    --------
    median

    """

    a = np.array(a, copy=False)
    a_median = np.median(a, axis=axis)

    #re-broadcast the output median array to subtract it
    if axis is not None:
        shape = list(a_median.shape)
        shape.append(1)
        a_median = a_median.reshape(shape)

    #calculated the median average deviation
    return np.median(np.abs(a - a_median), axis=axis)


def biweight_location(a, c=6.0, M=None):
    """
    Compute the biweight location for an array

    Returns the biweight location for the array elements.  The biweight
    is a robust statistic for determining the central location of a
    distribution.

    The biweight location is given by the follow equation::
    ..math::
        C_{bl}= M+\frac{\Sigma_{|u_i|<1} (x_i-M)(1-u_i^2)^2}
        {\Sigma_{|u_i|<1} (1-u_i^2)^2}
    where M is the sample mean or if run iterative the initial guess,
    and u_i is given by::
    ..math::a
        u_{i} = \frac{(x_i-M)}{cMAD}
    where MAD is the median absolute deviation.

    For more details, see Beers, Flynn, and Gebhardt, 1990, AJ, 100, 32B

    Parameters
    ----------
    a : array_like
        Input array or object that can be converted to an array.
    c : float
        Tuning constant for the biweight estimator.  Default value is 6.0.
    M : float, optional
        Initial gues for the biweight location.

    Returns
    -------
    biweight_location: float
        Returns the biweight location for the array elements.

    Examples
    --------

    This will generate random variates from a Gaussian distribution and return
    the median absolute deviation for that distribution::

        >>> from astropy.tools.alg import biweight_location
        >>> from numpy.random import randn
        >>> randvar = randn(10000)
        >>> cbl = biweight_location(randvar)

    See Also
    --------
    median absolute deviation, biweight_midvariance
    """

    a = np.array(a, copy=False)

    if M is None:
        M = np.median(a)

    #set up the difference
    d = a - M

    #set up the weighting
    u = d / c / median_absolute_deviation(a)
    

    #now remove the outlier points
    mask = np.abs(u) < 1
    u = (1 - u**2)**2

    return M+(d[mask]*u[mask]).sum()/u[mask].sum()


def biweight_midvariance(a, c=9.0, M=None):
    """
    Compute the biweight midvariance for an array

    Returns the biweight midvariance for the array elements.  The biweight
    midvariance is a robust statistic for determining the midvariance (ie. the
    standard deviation) of a distribution.

    The biweight location is given by the follow equation::
    ..math::
        C_{bl}= n^{1/2} \frac{[\Sigma_{|u_i|<1} (x_i-M)**2(1-u_i^2)^4]^{0.5}}
        {|\Sigma_{|u_i|<1} (1-u_i^2)(1-5u_i^2)|}
    where  u_i is given by::
    ..math::a
        u_{i} = \frac{(x_i-M)}{cMAD}
    where MAD is the median absolute deviation.  For the midvariance  parameter, c is
    typically uses a value of 9.0.

    For more details, see Beers, Flynn, and Gebhardt, 1990, AJ, 100, 32B

    Parameters
    ----------
    a : array_like
        Input array or object that can be converted to an array.
    c : float
        Tuning constant for the biweight estimator.  Default value is 9.0.
    M : float, optional
        Initial gues for the biweight location.

    Returns
    -------
    biweight_midvariance: float
        Returns the biweight midvariance for the array elements.

    Examples
    --------

    This will generate random variates from a Gaussian distribution and return
    the median absolute deviation for that distribution::

        >>> from astropy.tools.alg import biweight_midvariance
        >>> from numpy.random import randn
        >>> randvar = randn(10000)
        >>> scl = biweight_midvariance(randvar)

    See Also
    --------
    median absolute deviation, biweight_location
    """

    a = np.array(a, copy=False)
    n = len(a)

    if M is None:
        M = np.median(a)

    #set up the difference
    d = a - M

    #set up the weighting
    u = d / c / median_absolute_deviation(a)
   

    #now remove the outlier points
    mask = np.abs(u) < 1
    n=mask.sum()
    u = u**2
    return n**0.5 * (d[mask] * d[mask] * (1 - u[mask])**4).sum()**0.5 \
           / np.abs(((1 - u[mask]) * (1 - 5 * u[mask])).sum())
