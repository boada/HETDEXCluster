import numpy
from astLib import astCoords, astWCS
from astropy.io import fits
from itertools import count, izip
from scipy import optimize
from scipy.ndimage import zoom

def contour_levels(x, y=[], bins=10, levels=(0.68,0.95)):
    """
    Get the contour levels corresponding to a set of percentiles (given as
    fraction of 1) for a 2d histogram.

    Parameters
    ----------
        x : array of floats
            if y is given then x must be a 1d array. If y is not given then
            x should be a 2d array
        y : array of floats (optional)
            1d array with the same number of elements as x
        bins : argument of numpy.histogram2d
        levels : list of floats between 0 and 1
            the fractional percentiles of the data that should be above the
            returned values

    Returns
    -------
        level_values : list of floats, same length as *levels*
            The values of the histogram above which the fractional percentiles
            of the data given by *levels* are

    """
    if len(y) > 0:
        if len(x) != len(y):
            msg = 'Invalid input for arrays; must be either 1 2d array'
            msg += ' or 2 1d arrays'
            raise ValueError(msg)
    else:
        if len(numpy.array(x).shape) != 2:
            msg = 'Invalid input for arrays; must be either 1 2d array'
            msg += ' or 2 1d arrays'
            raise ValueError(msg)
    def findlevel(lo, hist, level):
        return 1.0 * hist[hist >= lo].sum()/hist.sum() - level
    if len(x) == len(y):
        hist, xedges, yedges = numpy.histogram2d(x, y, bins=bins)
        hist = numpy.transpose(hist)
        extent = (xedges[0], xedges[-1], yedges[0], yedges[-1])
    elif len(y) == 0:
        hist = numpy.array(x)
    level_values = [optimize.bisect(findlevel, hist.min(), hist.max(),
                                    args=(hist,l)) for l in levels]
    return level

def contours_external(ax, imgwcs, contourfile, levels, colors, lw=1):
    """
    Draw contours from contourfile in the frame of imgwcs.

    """
    contourwcs = astWCS.WCS(contourfile)
    contourdata = fits.getdata(contourfile)
    while len(contourdata.shape) > 2:
        contourdata = contourdata[0]
    # convert coords
    ny, nx = contourdata.shape
    xo, yo = contourwcs.pix2wcs(-1, -1)
    x1, y1 = contourwcs.pix2wcs(nx, ny)
    xo, yo = imgwcs.wcs2pix(xo, yo)
    x1, y1 = imgwcs.wcs2pix(x1, y1)
    contourdata = zoom(contourdata, 3, order=3)
    ax.contour(contourdata, levels, colors=colors, linewidths=lw,
               extent=(xo,x1,yo,y1))
    return

def corner(X, config=None, names='', labels=None, bins=20, bins1d=20,
           clevels=(0.68,0.95), contour_reference='samples',
           truths=None, truths_in_1d=False, truth_color='r',
           smooth=False, likelihood=None, likesmooth=1, colors='k', cmap=None,
           ls1d='-', ls2d='solid', style1d='curve', medians1d=True,
           percentiles1d=True, background=None, bweight=None, bcolor='r',
           alpha=0.5, limits=None,
           ticks=None, show_contour=True, top_labels=False, output='',
           verbose=False, **kwargs):
    """
    Do a corner plot (e.g., with the posterior parameters of an MCMC chain).
    Note that there may still be some issues with the tick labels.

    Parameters
    ----------
      X         : array-like
                  all posterior parameters. Can also be the outputs of
                  more than one chain, given as an array of arrays of models
                  (e.g., X = [[A1, B1, C1], [A2, B2, C2]])
      config    : str (optional - NOT YET IMPLEMENTED)
                  name of file containing any parameters whose default values
                  should be modified. Format of the file is two columns,
                  where the first is the parameter name as listed here,
                  and the second is the value for that parameter. If the
                  parameter takes a list of values they should be comma-
                  separated, and multiple entries semi-colon-separated.
                  For example, a file containing
                        bins             20
                        bins1d           50
                        colors       yellow
                        ticks     2,3,4;10,11,12;3.2,3.3,3.4
                  would only modify these parameters. Note that because of the
                  content of the 'ticks' parameter, the chain must be a
                  three-parameter model.
      names     : list of strings (optional)
                  Names for each of the chains. Will be used to show a legend
                  in the (empty) upper corner
      labels    : list of strings (optional)
                  names of the parameters
      bins      : int or array of ints (default 20)
                  Number of bins for the contours in the off-diagonal panels.
                  Should be one value per chain, one value per parameter,
                  or have shape (nchains,nparams)
      bins1d    : int or array of ints (default 20)
                  Number of bins for the histograms or curves in the diagonal
                  panels. Should be one value per chain, one value per
                  parameter, or have shape (nchains,nparams)
      clevels   : list of floats between 0 and 1 (default: (0.68,0.95))
                  percentiles at which to show contours
      contour_reference : {'samples', 'likelihood'} (default 'samples')
                  whether to draw contour on fractions of samples or
                  on likelihood levels. In the former case, *clevels*
                  must be floats between 0 and 1; in the latter, the
                  levels of the log-likelihood.
      truths    : {list of floats, 'medians', None} (default None)
                  reference values for each parameter, to be shown in
                  each panel
      smooth    : float (optional)
                  the width of the gaussian with which to smooth the
                  contours in the off-diagonal panels. If no value is given,
                  the contours are not smoothed.
      likelihood : array of floats (optional)
                  the likelihood surface, to be shown as a histogram in the
                  diagonals.
      likesmooth : int (default 1000)
                  the number of maxima to average over to show the
                  likelihood surface
      colors    : any argument taken by the *colors* argument of
                  pylab.contour(), or a tuple of them if more than one
                  model is to be plotted
      ls1d      : {'solid','dashed','dashdot','dotted'} (default 'solid')
                  linestyle for the diagonal plots, if style1d=='curve'.
                  Can specify more than one value as a list if more than one
                  model is being plotted.
      ls2d      : {'solid','dashed','dashdot','dotted'} (default 'solid')
                  linestyle for the contours. Can specify more than one value
                  as a list if more than one model is being plotted.
      style1d   : {'bar', 'step', 'stepfilled', 'curve'} (default 'curve')
                  if 'curve', plot the 1d posterior as a curve; else this
                  parameter is passed to the 'histtype' argument in
                  pyplot.hist()
      medians1d : bool (default True)
                  whether to show the medians in the diagonal panels as
                  vertical lines
      percentiles1d : bool (default True)
                  whether to show selected percentiles (see *clevels*) in the
                  diagonal panels as vertical lines
      background : {None, 'points', 'density', 'filled'} (default None)
                  If not None, then either points, a smoothed 2d histogram,
                  or filled contours are plotted beneath contours.
      bweight   : array-like, same length as e.g., A1
                  values to color-code background points
      bcolor    : color property, consistent with *background*
                  color of the points or filled contours, or colormap of the
                  2d density background. If truths are given they will be
                  shown in red and it is therefore recommended that the
                  colors be on a blue scale.
      alpha     : float between 0 and 1 (default 0.5)
                  transparency of the points if shown
      limits    : list of length-2 lists (optional)
                  a list of plot limits for each of the parameters.
      ticks     : list of lists (optional)
                  a list of tick positions for each parameter, to be printed
                  both in the x and y axes as appropriate.
      top_labels : boolean (default False)
                  whether to show axis and tick labels at the top of each
                  diagonal plot
      output    : string (optional)
                  filename to save the plot.
      verbose   : boolean
                  whether to print the marginalized values per variable
      kwargs    : keyword arguments to be passed to pylab.contour()


    Returns
    -------
      fig, axes_diagonal, axes_off : pylab figure and axes (diagonal and
                  off-diagonal) instances

    """
    import pylab
    from numpy import append, array, digitize, exp, histogram, histogram2d
    from numpy import linspace, median, percentile, sort, transpose
    from scipy.ndimage.filters import gaussian_filter
    if style1d == 'curve':
        from scipy import interpolate

    # not yet implemented
    options = _load_corner_config(config)
    # the depth of an array or list. Useful to assess the proper format of
    # arguments. Returns zero if scalar.
    depth = lambda L: (hasattr(L, '__iter__') and max(map(depth,L)) + 1) or 0

    nchains = (len(X) if depth(X) > 1 else 1)
    if nchains > 1:
        ndim = len(X[0])
        nsamples = len(X[0][0])
        if background == 'points':
            background = None
    else:
        ndim = len(X)
        nsamples = len(X[0])
        X = (X,)
    if nsamples == 0:
        msg = 'plottools.corner: received empty array.'
        msg += ' It is possible that you set the burn-in to be longer'
        msg += ' than the chain itself!'
        raise ValueError(msg)
    # check ticks
    if ticks is not None:
        if len(ticks) != ndim:
            print 'WARNING: number of tick lists does not match',
            print 'number of parameters'
            ticks = None
    # check limits
    if limits is not None:
        if len(limits) != ndim:
            print 'WARNING: number of limit lists does not match',
            print 'number of parameters'
            limits = None
    # check clevels - they should be fractions between 0 and 1.
    if 1 < max(clevels) <= 100:
        clevels = [cl/100. for cl in clevels]
    elif max(clevels) > 100:
        msg = 'ERROR: contour levels must be between 0 and 1 or between'
        msg += ' 0 and 100'
        print msg
        exit()
    # check truths
    if truths is not None:
        if len(truths) != ndim:
            truths = None
    # check likelihood
    if likelihood is not None:
        msg = 'WARNING: likelihood format not right - ignoring'
        lshape = likelihood.shape
        
        if len(lshape) == 1:
            likelihood = [likelihood]
        if lshape[0] != nchains or lshape[1] != nsamples \
            or len(lshape) != 2:
            print msg
            likelihood = None
    try:
        if len(smooth) != len(X[0]):
            print 'WARNING: number smoothing widths must be equal to',
            print 'number of parameters'
            smooth = [0 for i in X[0]]
    except TypeError:
        if smooth not in (False, None):
            smooth = [smooth for i in X[0]]
    # check the binning scheme.
    meta_bins = [bins, bins1d]
    for i, bname in enumerate(('bins','bins1d')):
        bi = meta_bins[i]
        # will fail if bi is a scalar
        try:
            bidepth = depth(bi)
        except TypeError:
            bidepth = 0
        # will be the same message in all cases below
        msg = 'ERROR: number of {0} must equal either number'.format(bname)
        msg += ' of chains or number of parameters, or have shape'
        msg += ' (nchains,nparams)'
        # this means binning will be the same for all chains
        ones = numpy.ones((nchains,ndim))
        # is it a scalar?
        if bidepth == 0:
            meta_bins[i] = bi * ones
        # or a 1d list?
        elif bidepth == 1:
            bi = numpy.array(bi)
            if len(bi) == ndim:
                meta_bins[i] = ones * bi
            elif len(bi) == nchains:
                meta_bins[i] = ones * bi[:,numpy.newaxis]
            else:
                print msg
                exit()
        elif (bidepth == 2 and nchains > 1 and \
              numpy.array(bi).shape != ones.shape) or \
             bidepth > 2:
            print msg
            exit()
    bins, bins1d = meta_bins
    # figure size
    if ndim > 3:
        figsize = 2 * ndim
    else:
        figsize= 3 * ndim
    axsize = 0.85 / ndim
    if len(X) == 1:
        if isinstance(colors, basestring):
            color1d = colors
        else:
            color1d = 'k'
    else:
        if len(colors) == len(X):
            color1d = colors
        # supports up to 12 names (plot would be way overcrowded!)
        else:
            color1d = ('g', 'orange', 'c', 'm', 'b', 'y',
                       'g', 'orange', 'c', 'm', 'b', 'y')
    if isinstance(ls1d, basestring):
        ls1d = [ls1d for i in X]
    if isinstance(ls2d, basestring):
        ls2d = [ls2d for i in X]
    # all set!
    axvls = ('--', ':', '-.')
    fig = pylab.figure(figsize=(figsize,figsize))
    # diagonals first
    plot_ranges = []
    axes_diagonal = []
    # for backward compatibility
    histtype = style1d.replace('hist', 'step')
    for i in xrange(ndim):
        ax = pylab.axes([0.1+axsize*i, 0.95-axsize*(i+1),
                         0.95*axsize, 0.95*axsize],
                        yticks=[])
        axes_diagonal.append(ax)
        if i < ndim-1:
            ax.set_xticklabels([])
        peak = 0
        edges = []
        for m, Xm in enumerate(X):
            edges.append([])
            if style1d == 'curve':
                ho, e = histogram(Xm[i], bins=bins1d[m][i], normed=True)
                xo = 0.5 * (e[1:] + e[:-1])
                xn = linspace(xo.min(), xo.max(), 500)
                n = interpolate.spline(xo, ho, xn)
                ax.plot(xn, n, ls=ls1d[m], color=color1d[m])
            else:
                n, e, patches = ax.hist(Xm[i], bins=bins1d[m][i],
                                        histtype=histtype,
                                        color=color1d[m], normed=True)
            edges[-1].append(e)
            if n.max() > peak:
                peak = n.max()
            area = n.sum()
            if medians1d:
                ax.axvline(median(Xm[i]), ls='-', color=color1d[m])
            if verbose:
                if len(names) == len(X):
                    print names[m]
                if labels is not None:
                    print '  %s' %(labels[i]),
                    if truths is None:
                        print ''
                    else:
                        print '({0})'.format(truths[i])
                    print ' ', median(Xm[i])
            for p, ls in izip(clevels, axvls):
                v = [percentile(Xm[i], 100*(1-p)/2.),
                        percentile(Xm[i], 100*(1+p)/2.)]
                if percentiles1d:
                    ax.axvline(v[0], ls=ls, color=color1d[m])
                    ax.axvline(v[1], ls=ls, color=color1d[m])
                if verbose:
                    print '    p%.1f  %.2f  %.2f' %(100*p, v[0], v[1])
        if likelihood is not None:
            for m, Xm, Lm, e in izip(count(), X, likelihood, edges):
                #print Lm.min(), Lm.max()
                binning = digitize(Xm[i], e[m])
                xo = 0.5 * (e[m][1:] + e[m][:-1])
                # there can be nan's because some bins have no data
                valid = array([(len(Lm[binning == ii]) > 0)
                               for ii in xrange(1, len(e[m]))])
                Lmbinned = [median(sort(Lm[binning == ii+1])[-likesmooth:])
                            for ii, L in enumerate(valid) if L]
                #Lmbinned = array(Lmbinned) + 100
                # normalized to the histogram area
                Lmbinned = exp(Lmbinned)
                Lmbinned -= Lmbinned.min()
                Lmbinned /= Lmbinned.sum() / area
                ax.plot(xo[valid], Lmbinned, '-',
                        color=truth_color, lw=3, zorder=-10)
        if truths_in_1d and truths is not None:
            ax.axvline(truths[i], ls='-', color=truth_color,
                       zorder=10)
        if i == ndim-1 and labels is not None:
            if len(labels) >= ndim:
                ax.set_xlabel(labels[i])
        # to avoid overcrowding tick labels
        if ticks is None:
            tickloc = pylab.MaxNLocator(4)
            ax.xaxis.set_major_locator(tickloc)
        else:
            ax.set_xticks(ticks[i])
        pylab.xticks(rotation=45)
        if limits is not None:
            ax.set_xlim(*limits[i])
        pylab.ylim(0, 1.1*peak)
        if i != ndim-1:
            ax.set_xticklabels([])
        if top_labels:
            topax = ax.twiny()
            topax.set_xlim(*ax.get_xlim())
            topax.xaxis.set_major_locator(tickloc)
            topax.set_xlabel(labels[i])
        plot_ranges.append(ax.get_xlim())
    # lower off-diagonals
    axes_off = []
    for i in xrange(1, ndim): # vertical axes
        for j in xrange(i): # horizontal axes
            ax = pylab.axes([0.1+axsize*j, 0.95-axsize*(i+1),
                             0.95*axsize, 0.95*axsize])
            axes_off.append(ax)
            extent = append(plot_ranges[j], plot_ranges[i])
            for m, Xm in enumerate(X):
                h, xe, ye = histogram2d(Xm[j], Xm[i], bins=bins[m][i])
                h = h.T
                extent = (xe[0], xe[-1], ye[0], ye[-1])
                if smooth not in (False, None):
                    h = gaussian_filter(h, (smooth[i],smooth[j]))
                levels = contour_levels(Xm[j], Xm[i], bins=bins[m][i],
                                        levels=clevels)
                if background == 'points':
                    if not (cmap is None or bweight is None):
                        ax.scatter(Xm[j], Xm[i], c=bweight, marker='.',
                                   s=4, lw=0, cmap=cmap, zorder=-10)
                    else:
                        ax.plot(Xm[j], Xm[i], ',',
                                color=bcolor, alpha=alpha, zorder=-10)
                elif background == 'density':
                    ax.imshow([Xm[i], Xm[j]], cmap=cm.Reds,
                               extent=extent)
                elif background == 'filled':
                    clvs = append(clevels, 1)
                    lvs = contour_levels(Xm[j], Xm[i], bins=bins[m][i],
                                         levels=clvs)
                    try:
                        if hasattr(bcolor[0], '__iter__'):
                            bcolor = [bc for bc in bcolor]
                    except TypeError:
                        pass
                    for l in xrange(len(levels), 0, -1):
                        if len(bcolor[l-1]) == 3:
                            bcolor[l-1] = [bcolor[l-1]]
                        ax.contourf(h, (lvs[l-1],lvs[l]),
                                    extent=extent, colors=bcolor[l-1])
                if show_contour:
                    ax.contour(h, levels, colors=color1d[m],
                               linestyles=ls2d[m], extent=extent,
                               zorder=10, **kwargs)
                if truths is not None:
                    #pylab.axvline(truths[j], ls='-', color=(0,0.5,1))
                    #pylab.axhline(truths[i], ls='-', color=(0,0.5,1))
                    ax.plot(truths[j], truths[i], '+',
                            color=truth_color, mew=4, ms=12, zorder=10)
            if labels is not None:
                if len(labels) == ndim:
                    if j == 0:
                        ax.set_ylabel(labels[i])
                    if i == ndim - 1:
                        ax.set_xlabel(labels[j])
            if j > 0:
                ax.set_yticklabels([])
            if i < ndim - 1:
                ax.set_xticklabels([])
            ax.set_xlim(*plot_ranges[j])
            ax.set_ylim(*plot_ranges[i])
            if ticks is not None:
                ax.set_xticks(ticks[j])
                ax.set_yticks(ticks[i])
            else:
                # to avoid overcrowding tick labels
                xloc = pylab.MaxNLocator(4)
                ax.xaxis.set_major_locator(xloc)
                yloc = pylab.MaxNLocator(4)
                ax.yaxis.set_major_locator(yloc)
            pylab.xticks(rotation=45)
    # dummy legend axes
    if len(X) > 1 and len(names) == len(X):
        lax = pylab.axes([0.1+axsize*(ndim-1), 0.95,
                          0.95*axsize, 0.95*axsize],
                         xticks=[], yticks=[])
        lax.set_frame_on(False)
        for c, model in izip(color1d, names):
            pylab.plot([], [], ls='-', lw=2, color=c, label=model)
        lg = pylab.legend(loc='center', ncol=1)
        lg.get_frame().set_alpha(0)
    if output:
        pylab.savefig(output, format=output[-3:])
        pylab.close()
    return fig, axes_diagonal, axes_off

def wcslabels(wcs, xlim, ylim, xsep='00:00:01', ysep='00:00:15'):
    """
    Get WCS ticklabels

    Parameters
    ----------
        wcs     : astWCS.WCS instance
                  the wcs of the image to be shown
        xlim    : sequence of length 2
                  the minimum and maximum values of the x axis
        ylim    : sequence of length 2
                  the minimum and maximum values of the y axis
        xsep    : string
                  separation of right ascension ticks in the x axis,
                  in colon-separated hms format
        xsep    : string
                  separation of declination ticks in the y axis, in
                  colon-separated dms format

    Returns
    -------
        [xticks, xticklabels] : lists containing the positions and labels
                  for right ascension hms labels
        [yticks, yticklabels] : lists containing the positions and labels
                  for declination dms labels

    """
    def roundout(label):
        if label[2] > 59:
            label[1] += 1
            label[2] -= 60
        if label[1] > 59:
            label[0] += 1
            label[1] -= 60
        return label
    left, right = xlim
    bottom, top = ylim
    wcslim = [wcs.pix2wcs(left, bottom), wcs.pix2wcs(right, top)]
    ralim, declim = numpy.transpose(wcslim)
    hmslim = [astCoords.decimal2hms(x, ':') for x in ralim]
    dmslim = [astCoords.decimal2dms(y, ':') for y in declim]
    if dmslim[0][0] == '-':
        sgn = -1
    else:
        sgn = 1
    # assumes that East is left, as usual
    xsep = numpy.array(xsep.split(':'), dtype=int)
    ysep = numpy.array(ysep.split(':'), dtype=int)
    xseconds = [float(h.split(':')[2]) for h in hmslim[::-1]]
    yseconds = [float(d.split(':')[2]) for d in dmslim]
    xlim = []
    ylim = []
    for i in xrange(2):
        xlim.append(hmslim[-i-1].split(':'))
        xlim[i] = [int(xlim[i][0]), int(xlim[i][1]), float(xlim[i][2])]
        ylim.append(dmslim[i].split(':'))
        ylim[i] = [int(ylim[i][0]), int(ylim[i][1]), float(ylim[i][2])]
    if dmslim[0][0] == '-':
        ylim = ylim[::-1]
    xticklabels = [numpy.array([int(x) for x in xlim[0]])]
    yticklabels = [numpy.array([int(y) for y in ylim[0]])]
    for i in xrange(3):
        if xsep[i] != 0:
            while xticklabels[0][i] % xsep[i] != 0:
                xticklabels[0][i] += 1
        if ysep[i] != 0:
            while yticklabels[0][i] % ysep[i] != 0:
                yticklabels[0][i] += 1
    xticklabels[0] = roundout(xticklabels[0])
    yticklabels[0] = roundout(yticklabels[0])
    while numpy.any(xticklabels[-1] + xsep < xlim[1]):
        xticklabels.append(xticklabels[-1] + xsep)
        xticklabels[-1] = roundout(xticklabels[-1])
    while numpy.any(yticklabels[-1] + ysep < ylim[1]):
        yticklabels.append(yticklabels[-1] + ysep)
        yticklabels[-1] = roundout(yticklabels[-1])
    for i in xrange(len(xticklabels)):
        xticklabels[i] = [('%2d' %x).replace(' ', '0')
                          for x in xticklabels[i]]
    for i in xrange(len(yticklabels)):
        yticklabels[i] = [('%2d' %y).replace(' ', '0')
                          for y in yticklabels[i]]
        if -10 < ylim[0][0] < 0:
            yticklabels[i][0] = '-0%d' %abs(int(yticklabels[i][0]))
    xticklabels = [':'.join(x) for x in xticklabels]
    yticklabels = [':'.join(y) for y in yticklabels]
    if dmslim[0][0] == '-':
        yticklabels = ['-{0}'.format(y) if y[0] != '-' else y
                       for y in yticklabels]
    x = [astCoords.hms2decimal(xtl, ':') for xtl in xticklabels]
    y = [astCoords.dms2decimal(ytl, ':') for ytl in yticklabels]
    xticks = [wcs.wcs2pix(i, min(declim))[0] for i in x]
    yticks = [wcs.wcs2pix(max(ralim), i)[1] for i in y]
    return [xticks, xticklabels], [yticks, yticklabels]

def _load_corner_config(config):
    """
    Not implemented!

    """
    options = {}
    # is there a configuration file at all!?
    if config is None:
        return options
    data = numpy.loadtxt(config, dtype=str, unpack=True)
    for key, value in izip(*data):
        values = value.split(';')
        ndim = len(values)
        values = [val.split(',') for val in values]
        for i in xrange(ndim):
            for j in xrange(len(values)):
                try:
                    values[i][j] = float(values[i][j])
                except ValueError:
                    pass
                try:
                    values[i][j] = int(values[i][j])
                except ValueError:
                    pass
        if ndim == 1:
            values = values[0]
        options[key] = values
    return options
