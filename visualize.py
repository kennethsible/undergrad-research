from matplotlib import pyplot as plt
from numpy import loadtxt, log, abs
import sys, warnings

data = loadtxt(sys.argv[1])
name = sys.argv[1].split('.')[0]
x, y = data[:, 0], data[:, 1]
plt.plot(x, y)
plt.savefig(name + '.png', dpi=250)

plt.figure()
warnings.filterwarnings("ignore")
plt.plot(log(x), log(abs(y)))
plt.xlim((6.5, 7))
plt.savefig(name + '_tail.png', dpi=250)
