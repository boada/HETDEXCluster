import pylab as pyl

d = pyl.genfromtxt('countsandMassErrors.txt', names=True)

f, (ax1, ax2) = pyl.subplots(1, 2)

bins = pyl.arange(11.5, 16, 0.5)
centers = (bins[:-1] + bins[1:]) / 2.
for i, (
        style, mode
) in enumerate(zip(['-', '--', '-.'], ['perfect', 'targeted', 'survey'])):
    N = pyl.array(d['counts'][i * 8:i * 8 + 8])
    if N[0] == 0.:
        massError = pyl.array(d['massErr'][i * 8 + 1:i * 8 + 8])
        N = pyl.array(d['counts'][i * 8 + 1:i * 8 + 8])
        ax1.plot(centers[1:],
                 1 / pyl.sqrt(N),
                 style,
                 label='Possion_' + mode,
                 marker='o')
        ax1.plot(centers[1:],
                 massError,
                 style,
                 label='MassErr_' + mode,
                 marker='o')
        ax2.plot(centers[1:],
                 massError**2 + 1 / N,
                 style,
                 label=mode,
                 marker='o')
    else:
        massError = pyl.array(d['massErr'][i * 8:i * 8 + 8])
        N = pyl.array(d['counts'][i * 8:i * 8 + 8])
        ax1.plot(centers,
                 1 / pyl.sqrt(N),
                 style,
                 label='Possion_' + mode,
                 marker='o')
        ax1.plot(centers,
                 massError,
                 style,
                 label='MassErr_' + mode,
                 marker='o')
        ax2.plot(centers, massError**2 + 1 / N, style, label=mode, marker='o')

ax1.semilogy()
ax2.semilogy()
