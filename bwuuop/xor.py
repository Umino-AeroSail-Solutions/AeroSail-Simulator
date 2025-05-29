import matplotlib.pyplot as plt
import numpy as np
import numpy.random as rnd
rnd.seed(20051101)
t = np.arange(0,2**8)
for i in range(2**8):
    r = rnd.randint(2**8)
    z = r^t
    plt.plot(t,z,)
plt.show()