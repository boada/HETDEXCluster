import os
import time

for i in range(100):
    print(i)
    os.system('bsub < jobscript')
    time.sleep(1)

