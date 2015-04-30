import os
import time

for i in range(17,101):
    print i
    os.system('bsub < jobscript')
    time.sleep(5)

