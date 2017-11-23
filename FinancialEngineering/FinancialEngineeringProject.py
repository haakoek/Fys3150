import numpy as np
import matplotlib.pyplot as plt
import sys

m0 = 1

m = np.linspace(0,100*m0,1000)

plt.plot(m,np.exp(-m))
plt.show()

