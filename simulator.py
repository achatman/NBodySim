"""
Author: Miles Lucas
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Particle(object):

    def __init__(self, r, v, m=1):
        self.r = np.array(r, dtype='float64')
        self.v = np.array(v, dtype='float64')
        self.m = m

    def __str__(self):
        return 'Pos: '+str(self.r)+', Vel: '+str(self.v)

    def potential(self, other):
        pot = -other.m * (self.r - other.r) / np.linalg.norm(self.r - other.r)**3
        return pot

    def update(self, a, dt):
        self.v += a*dt
        self.r += self.v*dt


class System(object):

    def __init__(self, bodies=[]):
        self.bodies = bodies

    def add(self, body):
        self.bodies.append(body)

    def update(self, dt):
        a = np.zeros(3)
        for b in self.bodies:
            for b2 in self.bodies:
                if (b != b2):
                    a += b.potential(b2)
            b.update(a, dt)

    def run(self, time, dt):
        fig = plt.figure()
        ax = Axes3D(fig)
        plt.style.use('classic')
        for i in range(time*60):
            # All for setting up axis
            ax.set_xlim(-10, 10)
            ax.set_xlabel('X')
            ax.set_ylim(-10, 10)
            ax.set_ylabel('Y')
            ax.set_zlim(-10, 10)
            ax.set_zlabel('Z')
            ax.tick_params(axis='both', bottom='off', top='off', right='off', left='off',
                labelbottom='off', labeltop='off', labelright='off', labelleft='off')

            #The actual money
            sys.update(dt)
            for b in sys.bodies:
                ax.scatter(b.r[0], b.r[1], b.r[2], 'b')
            #60 fps
            plt.pause(1/60)
            plt.cla()

#-------------------------------------------------------------------------------

# Driver Script
p = []
for i in range(10):
    pos = (np.random.rand(3)-0.5)*8
    vel = (np.random.rand(3)-0.5)*.5
    p.append(Particle(pos, vel))

sys = System(p)
sys.run(5, .1)
