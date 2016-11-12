import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import sys
import pylab

try:
	filename = sys.argv[1]
except:
	sys.exit('No argument. Enter valid file as argument')

try:
	f = open(filename, 'r')
except:
	sys.exit('No such file')

try:
	filename2 = sys.argv[2]
except:
	sys.exit('No argument. Enter valid file as argument')

try:
	f2 = open(filename2, 'r')
except:
	sys.exit('No such file')

f.readline()
sizeofgrid = int(f.readline())
f2.readline()
f2.readline()

intensity1 = [[0]*sizeofgrid for i in range(sizeofgrid)]
intensity2 = [[0]*sizeofgrid for i in range(sizeofgrid)]

for i in xrange(0, sizeofgrid):
	for j in xrange(0, sizeofgrid):
		intensity1[i][j] = float(f.readline())
		intensity2[i][j] = float(f2.readline())

#intensity1[16][16] = 0
#intensity2[16][16] = 0

plt.figure(1)

plt.subplot(211)
plt.title("Simulated picture of scattering through the interstellar medium cloud - Homogenous Medium")
plt.xlabel("X-pixel")
plt.ylabel("Y-pixel")
plt.imshow(intensity1, vmin=0, vmax=200, interpolation='nearest')
plt.grid(True)

plt.subplot(212)
plt.title("Simulated picture of scattering through the interstellar medium cloud - Nonhomogenous Medium")
plt.xlabel("X-pixel")
plt.ylabel("Y-pixel")
plt.imshow(intensity2, vmin=0, vmax=200, interpolation='nearest')
plt.grid(True)

plt.show()
