
"""
Author: Andrew Chatman
Last updated: 2017-04-17
"""

from Body import Body
from Body import DimensionException

"""
Represents a system of Body objects. A System
can simulate a dynamical N-Body graviational system.
"""
class System(object):
  
  """
  Initializes an empty system.
  """
  def __init__(self):
    self.numDim = -1
    self.bodies = []
    self.steps = 0
  
  """
  Adds a body to the system. Throws a DimensionException
  if the body does not exist in the same number of dimensions
  as the system.
  """
  def addBody(self, body):
    if(numDim == -1):
      numDim = len(body.pos)
    
    if(len(body.pos) != numDim):
      raise DimensionException("Bodies added must have the same dimensionality as the system.")
    
    bodies.append(body)
