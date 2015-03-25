#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pylab as pyl
import h5py as hdf

f = hdf.File('./stat_testing312170.hdf5','r')
dset = f[f.keys()[0]]
result = dset.value

# filter out bad values
x = pyl.where(result[:,2] != -1)
result = result[x[0]]
z = pyl.log10(result[:,1])
z = pyl.where(pyl.isnan(z) != True)
result = result[z[0]]

X = pyl.logspace(0.3, 3, 50)
Y = (result[:,2] - result[:,-3])/result[:,-3]

z = pyl.where(Y > 0)
M = Y[z[0]]

z = pyl.where(Y < 0)
W = Y[z[0]]

# Data to be represented
#X = pyl.linspace(0,100,20)
#M = (1+pyl.sin(X/20)+pyl.random.uniform(0.75,1.0,X.size)) * pyl.linspace(1,0.2,X.size)
#W = (1+pyl.sin(X/20)+pyl.random.uniform(0.75,1.0,X.size)) * pyl.linspace(1,0.2,X.size)

fig = pyl.figure(figsize=(8,6), dpi=72,facecolor="white")
axes = pyl.subplot(111, axisbelow=True)

index = pyl.digitize(result[:,1], bins) -1
avgs = [pyl.mean(y[index == k]) for k in range(len(bins))]
pyl.plot(bins-delta/2, avgs, lw=2, color=c, label=i)

pyl.bar(X, +M, width=5.25, facecolor='#9999ff', edgecolor='white')
for i,x in enumerate(X):
    x += 5.25/2
    y = M[i] + 0.1
    pyl.text(x,y,"%.1f" % y, color='#9999ff', size=9,
             horizontalalignment='center',  verticalalignment='bottom')

pyl.bar(X, -W, width=5.25, facecolor='#ff9999', edgecolor='white')
for i,x in enumerate(X):
    x += 5.25/2
    y = -W[i] - 0.1
    pyl.text(x,y,"%.1f" % (-y), color='#ff9999', size=9,
             horizontalalignment='center',  verticalalignment='top')

axes.set_xlim([0,110])
axes.set_xticks([])
axes.set_yticks([+1,-1])
axes.set_yticklabels(['MEN', 'WOMEN'])
labels = axes.get_yticklabels()
labels[0].set_rotation(90)
labels[0].set_color('#9999ff')
labels[1].set_rotation(90)
labels[1].set_color('#ff9999')

axes.spines['top'].set_color('none')
axes.spines['left'].set_color('none')
axes.spines['right'].set_color('none')
axes.spines['bottom'].set_color('none')
axes.yaxis.set_ticks_position('left')


pyl.text(108,2.00,"2012", color='black', size=36,
         horizontalalignment='right', verticalalignment='bottom')
pyl.text(108,1.95,"Those figures are\n totally random", color='.75', size=8,
         horizontalalignment='right', verticalalignment='top')

pyl.show()
