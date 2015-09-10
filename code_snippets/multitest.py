from multiprocessing import Pool
from time import sleep
from random import randint
import os


class AsyncFactory:
    def __init__(self, func, cb_func):
        self.func = func
        self.cb_func = cb_func
        self.pool = Pool(maxtasksperchild=2)

    def call(self,*args, **kwargs):
        self.pool.apply_async(self.func, args, kwargs, self.cb_func)

    def wait(self):
        self.pool.close()
        self.pool.join()



def square(pos, x):
    sleep_duration = randint(1,5)
    print "PID: %d \t Value: %d \t Sleep: %d" % (os.getpid(), x ,sleep_duration)
    sleep(sleep_duration)
    return pos, x*x

def cb_func((pos, x)):
    print pos, x
    results.append(x)


results = []

async_square = AsyncFactory(square, cb_func)

for i in range(25):
    async_square.call(i,i)

async_square.wait()
