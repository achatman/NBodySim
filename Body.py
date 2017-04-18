
"""
Author: Andrew Chatman
Last updated: 2017-04-17
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from random import randint
import math

#Constant. Please don't change.
G = 1; #6.67 * (10 ** -11);


"""
Dimension Exception is thrown when there is a 
mismatch in dimensionality, usually between two
or more bodies.
"""
class DimensionException(Exception):
  pass



"""
Represents a body with mass which can interact 
gravitationally with other bodies.
"""
class Body(object):
  
  """
  Default constructor. Sets position,
  velocity, and acceleration to 0 (3D).
  Mass is set to one (to avoid divide by 0).
  """
  def __init__(self):
    self.pos = [0,0,0]
    self.vel = [0,0,0]
    self.acc = [0,0,0]
    self.mass = 1
  
  """
  Mass constructor. Sets position,
  velocity, and acceleration to 0 (3D).
  Mass is given as the argument.
  """
  def __init__(self, given_mass, **kwargs):
    self.mass = given_mass
    self.force = []
    self.pos = []
    self.vel = []
    self.acc = []
    if "position" in kwargs:
      self.pos = kwargs["position"]
    if "velocity" in kwargs:
      self.vel = kwargs["velocity"]
    if "acceleration" in kwargs:
      self.acc = kwargs["acceleration"]
    
    ndim = max([len(self.pos),len(self.vel),len(self.acc)])
    if ndim == 0:
      self.pos = [0,0,0]
      self.vel = [0,0,0]
      self.acc = [0,0,0]
      self.force = [0,0,0]
    else:
      while len(self.pos) < ndim:
        self.pos.append(0)
      while len(self.vel) < ndim:
        self.vel.append(0)
      while len(self.acc) < ndim:
        self.acc.append(0)
      while len(self.force) < ndim:
        self.force.append(0)

  def __str__(self):
    return "Mass = %s, Position = %s"%(self.mass,self.pos)
  
  def __eq__(self,other):
    if self.mass != other.mass:
      return False
    if len(self.pos) != len(other.pos):
      return False
    for i in range(0,len(self.pos)):
      if self.pos[i] != other.pos[i]:
        return False
      if self.vel[i] != other.vel[i]:
        return False
      if self.acc[i] != other.acc[i]:
        return False
    return True
  
  def __ne__(self,other):
    return not self == other
  
  """
  Returns the magnitude of the distance this body is away from the origin.
  """
  def distanceFromOrigin(self):
    squares = 0
    for g in pos:
      squares += g**2
    return math.sqrt(squares)
  
  """
  Returns the magnitude of this body's velocity.
  """
  def getVelMag(self):
    sq = 0
    for g in self.vel:
      sq += g*g
    return math.sqrt(sq)
  
  """
  Returns the magnitude of this body's acceleration.
  """
  def getAccMag(self):
    sq = 0
    for g in self.acc:
      sq += g*g
    return math.sqrt(sq)
  
  """
  Returns the distance between this body and other.
  Both bodies should exist in the same dimensional
  space, else an error is thrown.
  """
  def getDistanceMag(self,other):
    if(len(self.pos) != len(other.pos)):
      raise DimensionException("The dimensionalities of two bodies must be equal to calculate a distance.")
    distSq = 0;
    for i in range(0,len(self.pos)):
      distSq += (self.pos[i] - other.pos[i])**2
    return math.sqrt(distSq)
  
  """
  Returns the vectorial distance between this body and other.
  Raises DimensionException if not of equal dimensionalities.
  """
  def getDistance(self,other):
    if(len(self.pos) != len(other.pos)):
      raise DimensionException("The dimensionalities of two bodies must be equal to calculate a distance.")
    out = []
    for i in range(0,len(self.pos)):
      out.append(other.pos[i] - self.pos[i]);
    return out;
  
  """
  Returns the vector of the gravitational
  force between this body and other.
  """
  def getGravForce(self, other):
    distVec = self.getDistance(other)
    distMag = self.getDistanceMag(other)
    force = (G * self.mass * other.mass) / (distMag ** 2)
    distUnitVec = []
    forceVec = []
    for i in range(0,len(distVec)):
      distUnitVec.append(distVec[i] / distMag)
      forceVec.append(force * distUnitVec[i])
    return forceVec
    
  """
  Returns the gravitational potential energy between
  this body and other.
  """
  def getGravPotentialEnergy(self, other):
    energy = (G * self.mass * other.mass) / self.getDistanceMag(other)
    return energy
  
  """
  Returns the gravitational potential due to other
  at the position of this body.
  """
  def getGravPotential(self, other):
    pot = (G * other.mass) / self.getDistanceMag(other)
    return pot
  
  """
  Returns the kinetic energy of this body.
  """
  def getKineticEnergy(self):
    energy = .5 * self.mass * (self.getVelMag()**2)
    return energy
  
  """
  Adds the given force to the current force on this body.
  """
  def addForce(self, f):
    if(len(f) != len(self.pos)):
      raise DimensionException("The dimensionality of the given force must be the same as the body on which it acts.")
    for i in range(0,len(f)):
      self.force[i] += f[i]
  
  """
  Kinematically adjust the body's position, velocity, and acceleration based on the 
  current force. "Consumes" the current force.
  """
  def step(self, dt):
    central = True
    for g in self.pos:
      if g != 0:
        central = False
    if central:
      return
    ndim = len(self.pos)
    for i in range(0,ndim):
      self.pos[i] += self.vel[i] * dt
      self.vel[i] += self.acc[i] * dt
      self.acc[i] += self.force[i] / self.mass
      self.force[i] = 0

"""
A model of a system of particles capable of simulating gravitational
interactions in time steps.
"""
class System(object):
  """
  Initializes an empty system of dimension d.
  """
  def __init__(self,d):
    self.ndim = d
    self.particles = []

  def addBody(self, mass, pos, vel = [], acc = []):
    if(len(pos) != self.ndim):
      raise DimensionException("The position of the particles must be of the same dimensionality as the system.")
    if len(vel) == 0:
      while len(vel) < len(pos):
        vel.append(0)
    if len(acc) == 0:
      while len(acc) < len(pos):
        acc.append(0)
    self.particles.append(Body(mass,position = pos, velocity  = vel, acceleration = acc))
  
  def __str__(self):
    out = ""
    for g in self.particles:
      out += str(g) + '\n'
    return out

  """
  Steps the entire system over time dt.
  Returns a list of position vectors for each particle.
  """
  def step(self, dt):
    for g in self.particles:
      force = []
      for h in range(0,self.ndim):
        force.append(0)
      for h in self.particles:
        if g != h:
          g.addForce( g.getGravForce(h) )
    
    for g in self.particles:
      g.step(dt);
      
    out = []
    for g in self.particles:
      out.append(np.asarray(g.pos))
    return np.asarray(out)


#EXECUTION SCRIPT

def getColor(num):
  chars = "0123456789abcdef"
  red = 200
  redness = 256 * num / red
  blueness = 256 - redness
  red1 = redness / 16
  red2 = redness % 16
  blue1 = blueness / 16
  blue2 = blueness % 16
  outstr = chars[red1] + chars[red2]    \
  + "00"                                \
  + chars[blue1] + chars[blue2]
  return outstr

#Make system
"""
N = int(input("Number of particles: "))
n_frames = int(input("Step Number: "))
sys = System(3)

for i in range(N):
  pos = [randint(-50,50),randint(-50,50),randint(-50,50)]
  mass = 1
  sys.addBody(mass,pos)
sys.addBody(100,[0,0,0])
"""
#Solar System
n_frames = 200
sys = System(3)
sys.addBody(333000,[0,0,0])                   #Sol
sys.addBody(.05,[.387,0,0],[0,2581,0])    #Mercury
sys.addBody(.8,[.72,0,0],[0,645,0])    #Venus
sys.addBody(1,[1,0,0],[0,577,0])         #Earth
sys.addBody(.1,[1.5,0,0],[0,471,0])    #Mars
'''
sys.addBody(318,[5.2,0,0],[0,253,0])  #Jupiter
sys.addBody(95,[9.5,0,0],[0,187,0])    #Saturn
sys.addBody(14.5,[19,0,0],[0,132,0])   #Uranus
sys.addBody(17.1,[30,0,0],[0,105,0])   #Neptune
'''



with open("system.log", mode = 'w') as w:
  w.write(str(sys))

rt = []
for t in range(n_frames):
  rt.append(sys.step(.000001))
rt = np.asarray(rt)

#Make plot
fig = plt.figure()
ax = fig.add_axes([0,0,1,1], projection = '3d')
plt.style.use('ggplot')

ax.set_xlim(-2,2)
ax.set_xlabel('X')
ax.set_ylim(-2,2)
ax.set_ylabel('Y')
ax.set_zlim(-2,2)
ax.set_zlabel('Z')
ax.tick_params(axis = 'both', bottom = 'off', top = 'off', right = 'off', left = 'off',
               labelbottom = 'off', labeltop = 'off', labelright = 'off', labelleft = 'off')

pts = sum([ax.plot([], [], [], 'bo') for n in range(len(sys.particles))], [])

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

anim = animation.FuncAnimation(fig, animate, init_func = init, frames = n_frames, interval = 20, blit = True)
plt.show()