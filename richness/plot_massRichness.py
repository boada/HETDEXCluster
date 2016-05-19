from __future__ import division
import h5py as hdf
from scipy import stats
import pylab as pyl
from astLib import astStats
from line import fit

def mklogMass(theta):
    ''' Makes Log10 Mass from a given richness following the simet2016
    mass-richness relationship.

    '''
    m0, alpha, lambda0 = theta
    return lambda x: m0 + alpha*pyl.log10(x/lambda0)

def mklambda(theta):
    ''' Makes Lambda given a Log10 Mass. If mass is not Log10 it won't work
    correctly.

    '''

    m0, alpha, lambda0 = theta
    return lambda y: lambda0*(10**y/10**m0)**(1/alpha)

def add_subplot_axes(ax,rect,axisbg='w'):
    ''' Add an inset axes to the first subplot. '''
    fig = pyl.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax

# normalization, power-law index, and lambda0 of the lambda-mass relation
truth = 14.19312, 1.31, 30

### Fake Data for Testing ###
#N = 100
#x_true = np.linspace(1,250, N)
#y_true = mklogMass(truth)(x_true)

#x_err, y_err = 0, 0.250
#x_obs = stats.norm(x_true, x_err).rvs(N)
#y_obs = stats.norm(y_true, y_err).rvs(N)
with hdf.File('./targetedRealistic_MLmasses_corrected.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    target_mlMasses = dset['HALOID', 'ML_pred_3d', 'ML_pred_3d_err', 'M200c']

# mask out the values with failed ML masses
mask = (target_mlMasses['ML_pred_3d'] != 0)
target_mlMasses = target_mlMasses[mask]

with hdf.File('./surveyCompleteRealistic_MLmasses_corrected.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    survey_mlMasses = dset['HALOID', 'ML_pred_3d', 'ML_pred_3d_err', 'M200c']

# mask out the values with failed ML masses
mask = (survey_mlMasses['ML_pred_3d'] != 0)
survey_mlMasses = survey_mlMasses[mask]

f = pyl.figure(1, figsize=(7*(pyl.sqrt(5.)-1.0)/2.0,7))
ax1s = pyl.subplot2grid((3,1), (2,0))
ax1 = pyl.subplot2grid((3,1), (0,0), rowspan=2)

scatter = 0.25

for m, c, style, zo in zip([target_mlMasses, survey_mlMasses], ['#7A68A6',
    '#188487'], ['-', '--'], [1,2]):
    bins = pyl.arange(10,150,20)

    # add the noise to the true masses-- 0.25 dex at the moment
    m_obs = stats.norm(m['M200c'], scatter).rvs(m.size)
    # use the noisy masses to calculate an observed lambda
    lam_obs = mklambda(truth)(m_obs)

    y_ = astStats.runningStatistic(lam_obs, m['ML_pred_3d'],
            pyl.percentile, binNumber=bins, q=[16, 50, 84])
    quants = pyl.array(y_[1])
    ax1.plot(y_[0],quants[:,1], style, c=c, zorder=zo)
    ax1.fill_between(y_[0], quants[:,2], quants[:,0], facecolor=c,
            alpha=0.4, edgecolor=c)

    # now we are going to do the scatter panel on the bottom
    bins = pyl.arange(10,150,10)
    x = lam_obs
    y = m['ML_pred_3d']
    index = pyl.digitize(x, bins)
    print [y[index==k].size for k in range(1,bins.size)]
    running = [pyl.std(y[index==k]) for k in range(1,bins.size)]
    running = pyl.array(running)
    tmp = pyl.where(~pyl.isnan(running))
    print(pyl.mean(running[tmp[0]]))
    centers = (bins[:-1] + bins[1:])/2.
    ax1s.plot(centers, running, style, drawstyle='steps', c=c)
    #ax1s.axhline(pyl.mean(running[tmp[0]]), zorder=0, lw=1)


### Add Legend ###
##################
line1 = pyl.Line2D([], [], ls='-', color='#7A68A6')
line2 = pyl.Line2D([], [], ls='--', color='#188487')
ax1.legend((line1, line2), ('Targeted', 'Survey'), loc=2)

ax1s.axhline(scatter, zorder=0)

ax1s.set_xlabel('Richness, $\lambda$')
ax1.set_ylabel('Log $M_{pred, corr}$')

ax1s.set_ylabel('$\sigma_{M|\lambda}$ (dex)')
ax1.set_xticklabels([])
ax1s.set_yticks([0.1, 0.2, 0.3, 0.4])
ax1s.set_ylim(0,0.5)
ax1.set_xticks(pyl.arange(10,200,40).tolist())
ax1s.set_xticks(pyl.arange(10,200,40).tolist())
ax1s.set_xlim(1,200)
ax1.set_xlim(1,200)

#####################
#### Now all of the stuff for the insert axes ####
#####################

target_rec = []
survey_rec = []
for scatter in pyl.arange(0.1, 0.5, 0.05):
    for m, s in zip([target_mlMasses, survey_mlMasses], [0,1]):
        # add the noise to the true masses
        m_pert = stats.norm(m['M200c'], scatter).rvs(m.size)
        # use the noisy masses to calculate an observed lambda
        lam = mklambda(truth)(m_pert)

        lam_bins = pyl.arange(10, 140, 10)
        idx = pyl.digitize(lam, lam_bins)
        scatters = pyl.zeros(lam_bins.size-1)
        for i in range(1, lam_bins.size):
            scatters[i-1] = pyl.std(m['ML_pred_3d'][idx==i])

        sc = pyl.mean(scatters)

        if s:
            survey_rec.append(sc)
        else:
            target_rec.append(sc)

# add the insert
rect = [.6, .15, .4, .4]
inset = add_subplot_axes(ax1, rect)
inset.plot([0,1], [0,1], c='k', zorder=0)
inset.plot(pyl.arange(0.1, 0.5, 0.05),target_rec, '-', c='#7a68a6')
inset.plot(pyl.arange(0.1, 0.5, 0.05),survey_rec, '--', c='#188487')
inset.set_xlabel('$\sigma_{true}$', fontsize=12)
inset.set_ylabel('$\sigma_{rec}$', fontsize=12)
inset.set_xlim(0,0.5)
inset.set_ylim(0,0.5)

inset.set_xticks([0.1, 0.2, 0.3, 0.4])
inset.set_yticks([0.1, 0.2, 0.3, 0.4])


pyl.show()





