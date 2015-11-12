import pylab as pyl
from matplotlib.ticker import NullFormatter
import h5py as hdf
from sklearn.cross_validation import train_test_split

# load the data
f = hdf.File('./result_FullKnowledge.hdf5', 'r')
dset = f[f.keys()[0]]
truth = dset.value
f.close()
f = hdf.File('./result_targetedIdeal.hdf5', 'r')
dset = f[f.keys()[0]]
target = dset.value
f.close()
f = hdf.File('./surveyComplete.hdf5', 'r')
dset = f[f.keys()[0]]
survey = dset.value

# only use the results we actually have results for
mask = target['NGAL'] >=5
maskedTruth = truth[mask]
maskedTarget = target[mask]
maskedSurvey = survey[mask]

Ngrid = 41
grid = pyl.linspace(12., 16., Ngrid + 1)

x = pyl.log10(maskedTarget['M200c'])
y = pyl.log10(maskedTarget['MASS'])

mask = (13.5 < y) & (y <16)
y = y[mask]
x = x[mask]

X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.2)
x, y = X_train, y_train

H, xbins, ybins = pyl.histogram2d(x, y, grid)
H = H.T
H /= pyl.sum(H)

#------------------------------------------------------------
# plot the result
fig = pyl.figure()

# define axes
ax_Pxy = pyl.axes((0.2, 0.34, 0.27, 0.52))
ax_Px = pyl.axes((0.2, 0.14, 0.27, 0.2))
ax_Py = pyl.axes((0.1, 0.34, 0.1, 0.52))
ax_cb = pyl.axes((0.48, 0.34, 0.01, 0.52))
ax_Px_y = [pyl.axes((0.65, 0.62, 0.32, 0.23)),
           pyl.axes((0.65, 0.38, 0.32, 0.23)),
           pyl.axes((0.65, 0.14, 0.32, 0.23))]

# set axis label formatters
ax_Px_y[0].xaxis.set_major_formatter(NullFormatter())
ax_Px_y[1].xaxis.set_major_formatter(NullFormatter())

ax_Pxy.xaxis.set_major_formatter(NullFormatter())
ax_Pxy.yaxis.set_major_formatter(NullFormatter())

ax_Px.yaxis.set_major_formatter(NullFormatter())
ax_Py.xaxis.set_major_formatter(NullFormatter())

# draw the joint probability
pyl.axes(ax_Pxy)
H *= 1000
pyl.imshow(H, interpolation='nearest', origin='lower', aspect='auto',
           extent=[12, 16, 12., 16], cmap=pyl.cm.binary)

cb = pyl.colorbar(cax=ax_cb)
cb.set_label('$p(M_{True}, M_{Rec})$')
pyl.text(0, 1.02, r'$\times 10^{-3}$',
         transform=ax_cb.transAxes)

# draw p(x) distribution
ax_Px.plot(xbins[1:], H.sum(0), '-k', drawstyle='steps')

# draw p(y) distribution
ax_Py.plot(H.sum(1), ybins[1:], '-k', drawstyle='steps')

# tweak axis labels
ax_Px.set_xticks([12,13,14,15,16])
ax_Py.set_yticks([12,13,14,15,16])

# define axis limits
ax_Pxy.set_xlim(12, 16)
ax_Pxy.set_ylim(12, 16)
ax_Px.set_xlim(12, 16)
ax_Py.set_ylim(12, 16)

# label axes
ax_Pxy.set_xlabel('$M_{True}$')
ax_Pxy.set_ylabel('$M_{Rec}$')
ax_Px.set_xlabel('$M_{True}$')
ax_Px.set_ylabel('$p(M_{True})$')
ax_Px.yaxis.set_label_position('right')
ax_Py.set_ylabel('$M_{Rec}$')
ax_Py.set_xlabel('$p(M_{Rec})$')
ax_Py.xaxis.set_label_position('top')

# draw marginal probabilities
iy = pyl.digitize(X_test[:3], ybins)

colors = 'rgc'
axis = ax_Pxy.axis()
for i in range(len(iy)):
    # overplot range on joint probability
    ax_Pxy.axhspan(ybins[iy[i]], ybins[iy[i]+1], alpha=0.25, fc=colors[i])
    Px_y = H[iy[i]] / H[iy[i]].sum()
    ax_Px_y[i].plot(xbins[1:], Px_y, drawstyle='steps', c=colors[i])
    ind = pyl.argsort(Px_y)
    ax_Px_y[i].axvline(xbins[ind[-1]], label='predicted')
    ax_Px_y[i].axvline(y_test[i], c='r', label='true')
    ax_Px_y[i].yaxis.set_major_formatter(NullFormatter())
    ax_Px_y[i].set_ylabel('$p(M_{True} | %.1f)$' % ybins[iy[i]])
ax_Pxy.axis(axis)

# make the predictions for the test set_xlabel
pred = []
iy = pyl.digitize(y_test, ybins)
for i in iy:
    if i > H.shape[0]:
        pass
    else:
        Px_y = H[i] / H[i].sum()
        norm = pyl.sum(Px_y * pyl.diff(xbins))
        centers = (xbins[:-1] + xbins[1:])/2.
        mean = pyl.sum(centers * Px_y * pyl.diff(xbins))/norm
            #pred.append(xbins[ind])
        pred.append(mean)

ax_Px_y[2].set_xlabel('$M_{True}$')
ax_Px_y[2].set_xticks([12,13,14,15,16])
ax_Px_y[2].legend()

ax_Pxy.set_title('Joint Probability')
ax_Px_y[0].set_title('Conditional Probability')

pyl.show()

