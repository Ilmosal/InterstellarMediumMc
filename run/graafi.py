import matplotlib.pyplot as plt
import numpy as np
import sys

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


packs = int(f.readline())
f2.readline()

theta = []
phi = []
intensity = []

theta2 = []
phi2 = []
intensity2 = []


for i in xrange(0, packs):
	phi.append(round(float(f.readline()), 4))
	theta.append(round(float(f.readline()), 4))
	intensity.append(round(float(f.readline()), 4))

	phi2.append(round(float(f2.readline()), 4))
	theta2.append(round(float(f2.readline()), 4))
	intensity2.append(round(float(f2.readline()), 4))



plt.figure(1)

plt.subplot(221)
plt.hist(phi, 100, weights = intensity)
plt.title("Angular Distribution of Scattered Photons from a Homogenous Medium")
plt.xlabel("Phi")
plt.ylabel("Amount of photons")

plt.subplot(222)
plt.hist(theta, 100, weights = intensity)
plt.title("Angular Distribution of Scattered Photons from a Homogenous Medium")
plt.xlabel("cos(Theta)")
plt.ylabel("Amount of photons")

plt.subplot(223)
plt.hist(phi2, 100, weights = intensity2)
plt.title("Angular Distribution of Scattered Photons from a Non-Homogenous Medium")
plt.xlabel("Phi")
plt.ylabel("Amount of photons")


plt.subplot(224)
plt.hist(theta2, 100, weights = intensity2)
plt.title("Angular Distribution of Scattered Photons from a Non-Homogenous Medium")
plt.xlabel("cos(Theta)")
plt.ylabel("Amount of photons")

plt.show()
