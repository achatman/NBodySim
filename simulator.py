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
        pot = other.m * (self.r - other.r) / np.linalg.norm(self.r - other.r)**3
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
        return np.asarray(r)

#-------------------------------------------------------------------------------
N = 2
n_frames = 500
# Driver Script
p = []
for i in range(N):
    pos = (np.random.rand(3)-0.5)*100
    # vel = (np.random.rand(3)-0.5)*.5
    vel = np.zeros(3)
    p.append(Particle(pos, vel))

sys = System(p)

center = Particle([0, 0, 0], [0, 0, 0], 100)
sys.add(center)

rt = []
for t in range(n_frames):
    rt.append(sys.unpack())
    sys.update(.1)
rt = np.asarray(rt)

#-------------------------------------------------------------------------------
# Animation

fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1], projection='3d')
plt.style.use('ggplot')
# All for setting up axis
ax.set_xlim(-100, 100)
ax.set_xlabel('X')
ax.set_ylim(-100, 100)
ax.set_ylabel('Y')
ax.set_zlim(-100, 100)
ax.set_zlabel('Z')
ax.tick_params(axis='both', bottom='off', top='off', right='off', left='off',
    labelbottom='off', labeltop='off', labelright='off', labelleft='off')

pts = sum([ax.plot([], [], [], 'bo') for n in range(len(sys.bodies))], [])

def init():
    for pt in pts:
        pt.set_data([], [])
        pt.set_3d_properties([])
    return pts

def animate(i):
    for r, pt in zip(rt[i], pts):
        pt.set_data(r[0], r[1])
        pt.set_3d_properties(r[2])

    fig.canvas.draw()
    return pts

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=n_frames, interval=20, blit=True)
plt.show()
