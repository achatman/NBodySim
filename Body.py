
"""
Author: Andrew Chatman
Last updated: 2017-04-17
"""


#Constant. Please don't change.
G = 6.67 * (10 ** -11);


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
  def __init__(self,given_mass):
    self.pos = [0,0,0]
    self.vel = [0,0,0]
    self.acc = [0,0,0]
    self.mass = given_mass
    
  """
  Full constructor. Accepts all instance
  variables as arguments. Position, velocity, and
  acceleration should be given as lists and of the
  same dimension. Coordinates that do not have enough
  dimensions are padded with 0.
  """
  def __init__(self, pos, vel, acc, mas):
    d = max([len(pos),len(vel),len(acc)])
    while len(pos) < d:
      pos.append(0)
    while len(vel) < d:
      vel.append(0)
    while len(acc) < d:
      acc.append(0)
    self.pos = pos
    self.vel = vel
    self.acc = acc
    self.mass = mas
  
  """
  Returns the magnitude of the distance this body is away from the origin.
  """
  def distanceFromOrigin(self):
    squares = 0
    for g in pos:
      squares += g**2
    return sqrt(squares)
  
  """
  Returns the magnitude of this body's velocity.
  """
  def getVelMag(self):
    sq = 0
    for g in self.vel:
      sq += g*g
    return sqrt(sq)
  
  """
  Returns the magnitude of this body's acceleration.
  """
  def getAccMag(self):
    sq = 0
    for g in self.acc:
      sq += g*g
    return sqrt(sq)
  
  """
  Returns the distance between this body and other.
  Both bodies should exist in the same dimensional
  space, else an error is thrown.
  """
  def getDistance(self,other):
    if(len(self.pos) != len(other.pos)):
      raise new DimensionException("The dimensionalities of two bodies must be equal to calculate a distance.")
    distSq = 0;
    for i in range(0,len(self.pos)):
      distSq += self.pos[i] * other.pos[i]
    return sqrt(distSq)
  
  """
  Returns the magnitude of the gravitational
  force between this body and other.
  """
  def getGravForce(self, other):
    force = (G * self.mass * other.mass) / (getDistance(self,other) ** 2)
    return force
    
  """
  Returns the gravitational potential energy between
  this body and other.
  """
  def getGravPotentialEnergy(self, other):
    energy = (G * self.mass * other.mass) / getDistance(self,other)
    return energy
  
  """
  Returns the gravitational potential due to other
  at the position of this body.
  """
  def getGravPotential(self, other):
    pot = (G * other.mass) / getDistance(self,other)
    return pot
  
  """
  Returns the kinetic energy of this body.
  """
  def getKineticEnergy(self):
    energy = .5 * self.mass * (getVelMag(self)**2)
    return energy
  