
from fenics import *
from mshr import *
import numpy as np

R=50

R0=10

rho0=0.1

tol=1E-14

#mesh=RectangleMesh(Point(0,0), Point(1, np.pi*2), 50, 50, 'right/left')
domain=Rectangle(Point(0,0), Point(1, np.pi*2))

mesh=generate_mesh(domain, 32)
