import pylab as pyl
import h5py as hdf
from astLib import astCalc as aca
from scipy.optimize import curve_fit

def correctFit(x, a):
    ''' This is the general form of the VD-mass relation that we are going to
        use. This function is called by the curve fitter to correct the A1D
        parameter for our specific data case.

    '''

    return a * (x/1e15)**(1/3.)


def wrapper_curveFit(x, y, a):
    ''' This is a simple wrapper which lets us figure out the bootstrap CI
        for the parameter given from correctFit.

    '''

    opt_parms, parm_cov = curve_fit(correctFit, x, y, p0=[a])
    return opt_parms[0]

def mass(vd, z, a):
    ''' Calculates the mass.

    '''

    return (vd/a)**(1/0.3361) * 1e15/(aca.H0*aca.Ez(z)/100)



f = hdf.File('./stat_testing312170.hdf5','r')
dset = f[f.keys()[0]]
result = dset.value

# filter out bad values
x = pyl.where(result[:,2] != -1)
result = result[x[0]]
result[:,1] = pyl.log10(result[:,1])

bins = pyl.linspace(13, 14.5, 20)
#bins = pyl.logspace(13, 14.5, 10)
delta = bins[1] - bins[2]


result[:,-2] = pyl.log10(result[:,-2])

f, ax = pyl.subplots(2,3)

ax = pyl.ravel(ax)
for j, i in enumerate(range(2, 13, 2)):


    slope = wrapper_curveFit(10**(result[:,-2]), result[:,i], 1082)

    print slope
    if slope == pyl.NaN:
        slope = 1082
    print slope

    slope = 1082.9
    m = [mass(vd, z, slope) for vd, z in zip(result[:,i], result[:,-1])]
    y = (m - 10**(result[:,-2]))/10**(result[:,-2])
    #y = (m - result[:,-2])/(result[:,-2])



    #ax[j].text(13.5, 0.25, str(slope), color='k', size=30)
    #y = (result[:,i] - result[:,-3])/result[:,-3]

    try:
        z = pyl.where(y > 0)
        index = pyl.digitize(result[:,-2][z[0]], bins) -1
        avgs = [pyl.median(y[z[0]][index == k]) for k in range(len(bins))]

        ax[j].bar(bins-delta/2, avgs, width=+delta, facecolor='#AAF0D1',
                edgecolor='white')

        for i,x in enumerate(bins):
            x += -1*delta/4
            Y = avgs[i] + 0.001
            #print Y
            ax[j].text(x,Y,"%2.f" % (Y*100), color='k', size=9,
                horizontalalignment='right', verticalalignment='bottom')
    except:
        pass

    try:
        z = pyl.where(y < 0)
        index = pyl.digitize(result[:,-2][z[0]], bins) -1
        avgs = [pyl.median(y[z[0]][index == k]) for k in range(len(bins))]
        ax[j].bar(bins-delta/2, avgs, width=+delta, facecolor='#AAF0D1',
                edgecolor='white')

        for i,x in enumerate(bins):
            x += -1*delta/4
            #print x
            Y = avgs[i] - 0.002
            ax[j].text(x,Y,"%2.f" % (-Y*100), color='k', size=9,
                horizontalalignment='right', verticalalignment='top')
    except:
        pass

    try:
        index = pyl.digitize(result[:,-2], bins) -1
        avgs = [pyl.median(y[index == k]) for k in range(len(bins))]
        ax[j].bar(bins-delta/4, avgs, width=+delta/2, facecolor='#9999ff',
            edgecolor='white')
    except:
        pass

    #ax[j].set_xlim(0.75, 3.25)
    ax[j].set_ylim(-0.5, 0.5)
    ax[j].set_xticks([])
    ax[j].set_yticks([])

    ax[j].spines['top'].set_color('none')
    ax[j].spines['left'].set_color('none')
    ax[j].spines['bottom'].set_color('none')
    ax[j].spines['right'].set_color('none')

'''
# Single plot annotations
arrowprops = dict(arrowstyle="-",
        connectionstyle="angle, angleA=90, angleB=180, rad=0")
ax[0].annotate('IQR', xy=(bins[-3]-delta, 0.15), xytext=(3.1, 0.10),
textcoords='data', arrowprops=arrowprops, rotation='vertical')

arrowprops = dict(arrowstyle="-",
        connectionstyle="angle, angleA=90, angleB=180, rad=0")
ax[0].annotate('Median', xy=(bins[-2]-delta/1.3, -0.025), xytext=(3.1, -0.10),
textcoords='data', arrowprops=arrowprops, rotation='vertical')

arrowprops = dict(arrowstyle="-",
        connectionstyle="angle, angleA=180, angleB=-90, rad=0")
ax[0].annotate('5\nMembers', xy=(bins[0], -0.16), xytext=(1.1, -0.21),
textcoords='data', arrowprops=arrowprops)

arrowprops = dict(arrowstyle="-",
        connectionstyle="angle, angleA=90, angleB=180, rad=0")
ax[0].annotate('1000\nMembers', xy=(bins[-1], -0.16), xytext=(2.75, -0.21),
textcoords='data', arrowprops=arrowprops, ha='right')
'''
pyl.show()
