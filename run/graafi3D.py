from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import sys

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

try:
	filename = sys.argv[1]
except:
	sys.exit('No argument. Enter valid file as argument')

try:
	f = open(filename, 'r')
except:
	sys.exit('No such file')

packs = int(f.readline())

theta = []
phi = []
intensity = []

for i in xrange(0, packs):
	phi.append(round(float(f.readline()), 4))
	theta.append(round(float(f.readline()), 4))
	intensity.append(round(float(f.readline()), 4))

hist, xedges, yedges = np.histogram2d(phi, theta, bins=(15,15))
xpos, ypos = np.meshgrid(xedges[:-1]+xedges[1:], yedges[:-1]+yedges[1:])

xpos = xpos.flatten()/2.
ypos = ypos.flatten()/2.
zpos = np.zeros_like (xpos)

dx = xedges [1] - xedges [0]
dy = yedges [1] - yedges [0]
dz = hist.flatten()

ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b', zsort='average')
plt.xlabel ("X")
plt.ylabel ("Y")

plt.show()

