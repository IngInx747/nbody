#coding=utf-8

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
from optparse import OptionParser

def readCommand(argv):
    usage = ""
    parser = OptionParser()
    parser.add_option("-f", "--input", dest="input", action="store")
    parser.add_option("-s", "--output", dest="output", action="store")
    parser.add_option("-n", "--particle", dest="particle", type="int", action="store")
    parser.add_option("-i", "--iterations", dest="iterations", type="int", action="store")
    options, junks = parser.parse_args(argv)
    if len(junks) != 0:
        raise Exception("Cannot understand: " + str(junks))
    args = dict()
    args["ifilename"] = options.input
    args["ofilename"] = options.output
    args["nBodies"] = options.particle
    args["nIters"] = options.iterations
    return args

def get_data(ifilename, nParticles, nIterations):

    '''
    Read binary file and turn raw data into format: Matrix[T][N][3]
    '''
    
    # get double
    data = np.fromfile(ifilename)
    # get float
    #data = np.fromfile(ifilename, np.float32)

    N = nParticles
    T = nIterations
    data = [[data[t*3*N + x*3 : t*3*N + (x+1)*3] for x in range(N)] for t in range(T)]

    return data

# Function that will be called for each frame of the animation
def update(t):
    update.t += 1
    if (update.t * 100) % nIters == 0:
        print(".", end='', file=sys.stdout, flush=True)
    # Reads a set of bodies into an array
    arr = np.empty(dtype=float, shape=(nBodies,3))
    for i in range(nBodies):
        arr[i,:] = data[t][i][:]

    points.set_xdata(arr[:,0])
    points.set_ydata(arr[:,1])
    points.set_3d_properties(arr[:,2]) # No set_zdata, se we use this

    return points,

# Read commands from console
args = readCommand(sys.argv[1:])

# Input file name
fname = args["ifilename"]
oname = 'nbody_simulation.mp4'
if (args["ofilename"]):
    oname = args["ofilename"]
# Reads header info (nBodies, nIters)
nBodies = args["nBodies"]
nIters = args["nIters"]

print("loading input file ", fname," ...")

# Opens input file
data = get_data(fname, nBodies, nIters)

# Adjusts marker size based on number of bodies in problem
marker = 1.0
if (nBodies > 1000):
    marker = 0.5
if (nBodies > 10000 ):
    marker = 0.1
if( nBodies > 100000 ):
    marker = 0.02

# Allocations array to hold a timestep
arr = np.empty(dtype=float, shape=(nBodies,3))

# Reads initial conditions
for i in range(nBodies):
    arr[i,:] = data[0][i][:]

# Create a 3D plot and initialize it with initial particle states
fig, ax = plt.subplots()
ax = p3.Axes3D(fig)

# Build Plot
points, = ax.plot3D(arr[:,0], arr[:,1], arr[:,2], 'wo', markersize=marker)
ax.set_ylim(-2.0, 2.0)
ax.set_xlim(-2.0, 2.0)
ax.set_zlim3d(-2.0, 2.0)
ax.set_facecolor('xkcd:black')
plt.axis('off')
ax.grid(False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])

update.t = -1

# Generate the animation
ani = animation.FuncAnimation(fig, update, nIters-2)

# Save .mp4 of the animation
# Bitrate may need to be increased for higher particle counts
ani.save(oname, fps=60, bitrate=500000, extra_args=["-s", "1280x720"])
#plt.show()

print("Amination completed")
