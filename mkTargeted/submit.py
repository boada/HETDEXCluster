import os
import time
#from data_handler import mkHalo
#import numpy as np

#halo = mkHalo()
#mask = (halo['m200c']/0.72 >= 1e13) & (halo['upid'] == -1)
#maskedHalo = halo[mask]
#hids, uniqueIdx = np.unique(maskedHalo['id'], return_index=True)

#interval = hids.size//10
interval = 40528//10

print interval

jobscript = \
        '''
## send stderr and stdout to the same file
#BSUB -o output.%J

## login shell to avoid copying env from login session
## also helps the module function work in batch jobs
#BSUB -L /bin/bash

## 30 minutes of walltime ([HH:]MM)
#BSUB -W 02:00

## numprocs
#BSUB -n 20

## 20 cores/node
#BSUB -R 'span[ptile=20]'

source /home/boada/.bashrc
'''

for i in range(10):
    with open('jobscript'+str(i), 'w') as f:
        f.write('## job name\n')
        f.write('#BSUB -J targetedRealistic'+str(i)+'\n')
        f.write(jobscript)
        if i == 9:
            f.write('time py targetedRealistic_async.py '+str(i*interval)+ \
                    ' '+str(i))
        else:
            f.write('time py targetedRealistic_async.py '+str(i*interval)+ \
                    ' '+str(i*interval+interval))

    print i
    os.system('bsub < jobscript'+str(i))
    time.sleep(2)

