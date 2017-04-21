'''
Author: Andrew Chatman
'''

import numpy as np

#Constants
G = 1


class Particle(Object):
  def __init__(self, mass = 1,
                     pos = [0,0],
                     vel = [0,0],
                     acc = [0,0]):
    self.m = mass
    self.r = np.array(pos, dtype='float64')
    self.v = np.array(vel, dtype='float64')
    self.a = np.array(acc, dtype='float64')
    self.F = np.array([0,0], dtype='float64')
  
  def __str__(self):
    return "m = %s, r = %s"%(self.m,self.r);
  
  def rmag(self):
    return np.linalg.norm(self.r)
  def vmag(self):
    return np.linalg.norm(self.v)
  def amag(self):
    return np.linalg.norm(self.a)
  
  #Gravitational force between self and other
  def gravForce(self, other):
    #vec[F] = (r/|r|) * (GMm)/(|r|^2)
    #vec[F] = (GMmr)/(|r|^3)
    #where r = r2 - r1
    
    F = ( (G * self.m * other.m *
             (other.r - self.r)) /
              (np.linalg.norm(other.r - self.r)**3) )
    
    return F;
  #Gravitational potential energy between self and other
  def gravPotentialEnergy(self, other):
    #U = GMm/|r|
    #r = r2 - r1
    
    U = G * self.m * other.m / np.linalg.norm(other.r - self.r)
    
    return U;
  
  #Graviational field due to self at point P
  def gravField(self, P):
    #g = (r/|r|) * (GM/|r|^2)
    #g = GMr/|r|^3
    #r = r2 - r1
    
    g = G * self.m * (P - self.r) / np.linalg.norm(P = self.r)**3
    
    return g;
  
  #Gravitational Potential due to self at point P
  def gravPotential(self, P):
    #V = GM/|r|
    #r = r2 - r1
    
    V = G * self.m / np.linalg.norm(P - self.r)
    
    return V;
  
  #Kinetic energy of self
  def kineticEnergy(self):
    #T = .5mv^2
    
    T = .5 * self.m * np.linalg.norm(self.v)
    
    return T;
  
  #Linear momentum of self
  def momentum(self):
    #p = mv
    
    p = self.m * self.v
    
    return p;
  
  #Angular momentum of self about point P
  def angularMomentum(self, P):
    #L = r x p
    #r = r_P - self.r
    
    L = np.cross(P - self.r, self.momentum)
    
    return L;

  #Applies current force to particle and updates the state variables kinematically
  def step(self, dt):
    self.r += self.v * dt + .5 * self.a * dt**2
    self.v += self.a * dt
    self.a += self.F / self.m
    self.F = np.array([0,0], dtype='float64')


class System(object):
  def __init__(self):
    self.particles = []
    self.time = 0
  
  def addParticle(self, mass = 1,
                        pos = [0,0],
                        vel = [0,0],
                        acc = [0,0]):
    self.particles.append(mass,pos,vel,acc)
  
  #Net gravitational force on the ith element of self.particles
  def netGravForce(self, i):
    netF = np.array([0,0], dtype='float64')
    this = self.particles[i]
    for g in self.particles:
      if g != this:
        netF += this.gravForce(g)
    
    return netF;
  
  #Net gravitational potential energy on the ith element of self.particles
  def netGravPotEnergy(self, i):
    netU = 0
    this = self.particles[i]
    for g in self.particles:
      if g != this:
        netU += this.gravPotentialEnergy(g)
    
    return netU;
  
  #Net graviational field at P
  def netGravField(self, P):
    netg = np.array([0,0], dtype='float64')
    for g in self.particles:
      netg += g.gravField(P)
    
    return netg;
  
  #Net graviational potential at P
  def netGravPot(self, P):
    netV = 0
    for g in self.particles:
      netV += g.gravPotential(P)
    
    return netV;
  
  #Net momentum of the system
  def netMomentum(self):
    netp = np.array([0,0], dtype='float64')
    for g in self.particles:
      netp += g.momentum();
    
    return netp;
  
  #Net angular momentum of the system about point P
  def netAngularMomentum(self, P):
    netL = np.array([0,0,0], dtype='float64')
    for g in self.particles:
      netL += g.angularMomentum(P)
    
    return netL;
  
  #Total kinetic energy of the system
  def totalKinetic(self):
    T = 0
    for g in self.particles:
      T += g.kineticEnergy()
      
    return T;
  
  #Total potential energy of the system
  def totalPotentialEnergy(self):
    U = 0
    for i in range(0, len(self.particles)):
      U += self.netGravPotEnergy(i)
    
    return U;
  
  #Lagrangian of the system
  def lagrangian(self):
    return self.totalKinetic() - self.totalPotentialEnergy();
  
  #Hamiltonian of the system
  def hamiltonian(self):
    return self.totalKinetic() + self.totalPotentialEnergy();
  
  def step(self, dt):
    for i, g in enumerate(self.particles):
      g.F = self.netGravForce(i)
    
    for g in self.particles:
      g.step()
    
    self.time += dt
      