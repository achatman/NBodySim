"""
Author: Miles Lucas
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

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

    def unpack(self):
        r = []
        for b in self.bodies:
            r.append(b.r)
        return np.array(r)

#-------------------------------------------------------------------------------
N = 100
# Driver Script
p = []
for i in range(N):
    pos = (np.random.rand(3)-0.5)*10
    # vel = (np.random.rand(3)-0.5)*.5
    vel = np.zeros(3)
    p.append(Particle(pos, vel))

sys = System(p)


#-------------------------------------------------------------------------------
# Animation

fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1], projection='3d')
plt.style.use('classic')
# All for setting up axis
ax.set_xlim(-10, 10)
ax.set_xlabel('X')
ax.set_ylim(-10, 10)
ax.set_ylabel('Y')
ax.set_zlim(-10, 10)
ax.set_zlabel('Z')
ax.tick_params(axis='both', bottom='off', top='off', right='off', left='off',
    labelbottom='off', labeltop='off', labelright='off', labelleft='off')

pts = sum([ax.plot([], [], [], 'bo') for n in range(N)], [])

def init():
    for pt in pts:
        pt.set_data([], [])
        pt.set_3d_properties([])
    return pts

def animate(i):
    sys.update(.1)
    for r, pt in zip(sys.unpack(), pts):
        pt.set_data(r[0], r[1])
        pt.set_3d_properties(r[2])

    fig.canvas.draw()
    return pts

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=100, interval=20, blit=True)
plt.show()
