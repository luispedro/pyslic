import numpy as np
from numpy import *

__all__ = ['nhistogram']

# This is necessary to work with multiple numpy installations.
# It uses the new behaviour of numpy.
def nhistogram(a, bins=10, range=None, normed=False, weights=None):
    """
    Compute the histogram of a set of data.

    Parameters
    ----------
    a : array_like
        Input data.
    bins : int or sequence of scalars, optional
        If `bins` is an int, it defines the number of equal-width
        bins in the given range (10, by default). If `bins` is a sequence,
        it defines the bin edges, including the rightmost edge, allowing
        for non-uniform bin widths.
    range : (float, float), optional
        The lower and upper range of the bins.  If not provided, range
        is simply ``(a.min(), a.max())``.  Values outside the range are
        ignored. Note that with `new` set to False, values below
        the range are ignored, while those above the range are tallied
        in the rightmost bin.
    normed : bool, optional
        If False, the result will contain the number of samples
        in each bin.  If True, the result is the value of the
        probability *density* function at the bin, normalized such that
        the *integral* over the range is 1. Note that the sum of the
        histogram values will often not be equal to 1; it is not a
        probability *mass* function.
    weights : array_like, optional
        An array of weights, of the same shape as `a`.  Each value in `a`
        only contributes its associated weight towards the bin count
        (instead of 1).  If `normed` is True, the weights are normalized,
        so that the integral of the density over the range remains 1.
        The `weights` keyword is only available with `new` set to True.

    Returns
    -------
    hist : array
        The values of the histogram. See `normed` and `weights` for a
        description of the possible semantics.
    bin_edges : array of dtype float
        Return the bin edges ``(length(hist)+1)``.


    See Also
    --------
    histogramdd

    Notes
    -----
    All but the last (righthand-most) bin is half-open.  In other words, if
    `bins` is::

      [1, 2, 3, 4]

    then the first bin is ``[1, 2)`` (including 1, but excluding 2) and the
    second ``[2, 3)``.  The last bin, however, is ``[3, 4]``, which *includes*
    4.

    Examples
    --------
    >>> np.histogram([1,2,1], bins=[0,1,2,3])
    (array([0, 2, 1]), array([0, 1, 2, 3]))

    """
    a = asarray(a)
    if weights is not None:
        weights = asarray(weights)
        if np.any(weights.shape != a.shape):
            raise ValueError, 'weights should have the same shape as a.'
        weights = weights.ravel()
    a =  a.ravel()

    if (range is not None):
        mn, mx = range
        if (mn > mx):
            raise AttributeError, \
                'max must be larger than min in range parameter.'

    if not iterable(bins):
        if range is None:
            range = (a.min(), a.max())
        mn, mx = [mi+0.0 for mi in range]
        if mn == mx:
            mn -= 0.5
            mx += 0.5
        bins = linspace(mn, mx, bins+1, endpoint=True)
    else:
        bins = asarray(bins)
        if (np.diff(bins) < 0).any():
            raise AttributeError, 'bins must increase monotonically.'

    # Histogram is an integer or a float array depending on the weights.
    if weights is None:
        ntype = int
    else:
        ntype = weights.dtype
    n = np.zeros(bins.shape, ntype)

    block = 65536
    if weights is None:
        for i in arange(0, len(a), block):
            sa = sort(a[i:i+block])
            n += np.r_[sa.searchsorted(bins[:-1], 'left'), \
                sa.searchsorted(bins[-1], 'right')]
    else:
        zero = array(0, dtype=ntype)
        for i in arange(0, len(a), block):
            tmp_a = a[i:i+block]
            tmp_w = weights[i:i+block]
            sorting_index = np.argsort(tmp_a)
            sa = tmp_a[sorting_index]
            sw = tmp_w[sorting_index]
            cw = np.concatenate(([zero,], sw.cumsum()))
            bin_index = np.r_[sa.searchsorted(bins[:-1], 'left'), \
                sa.searchsorted(bins[-1], 'right')]
            n += cw[bin_index]

    n = np.diff(n)

    if not normed:
        return n, bins
    db = array(np.diff(bins), float)
    return n/(n*db).sum(), bins
