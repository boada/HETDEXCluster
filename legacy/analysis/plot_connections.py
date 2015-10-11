from itertools import cycle
clist = rcParams['axes.color_cycle']

scatter(cluster[:,0],cluster[:,1])

for inds,c in zip(indices[[2,6,10]], cycle(clist)):
    for i in range(1,inds.shape[0]-1):
        plot([cluster[:,0][inds][0], cluster[:,0][inds][i]],
                [cluster[:,1][inds][0], cluster[:,1][inds][i]], c=c)

