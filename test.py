from matplotlib import pyplot as mp
import numpy as np

def gaussian(x, mu, sig):
    return np.exp(-(x - mu)**2. / (2 * sig**2.))

#for mu, sig in [(-1, 1), (0, 2), (2, 3)]:
mu=100.
sig=20.
mp.plot(gaussian(np.arange(200), mu, sig))

mp.show()