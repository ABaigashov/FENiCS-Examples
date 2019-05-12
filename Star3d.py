
from fenics import *
from mshr import *
import numpy as np

R=50

R0=10

rho0=0.1

domain=Sphere(Point(0,0,0),R)

mesh=generate_mesh(domain, 32)

plot(mesh)

V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition
u_D = Expression('(1 - x[0]*x[0] + x[1]*x[1] + x[2]*x[2])*np.heaviside(1-x[0]*x[0]-x[1]*x[1]-x[2]*x[2])', degree=2)

def boundary(x, on_boundary):
	return on_boundary

bc = DirichletBC(V, Constant(0), boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression('(1 - x[0]*x[0] + x[1]*x[1] + x[2]*x[2])*np.heaviside(1-x[0]*x[0]-x[1]*x[1]-x[2]*x[2],1)', degree=3)
a = dot(grad(u), grad(v))*dx
L = f*v*dx

# Compute solution
u = Function(V)
solve(a == L, u, bc)

# Save to file and plot
File("classic.pvd") << u

'''
V = FunctionSpace(mesh, 'P', 2)

bc = DirichletBC(V, Constant(0), 'on_boundary')


rho=rho0*(1-x[0]**2-x[1]**2-x[2]**2)*np.heaviside(1-x[0]**2-x[1]**2-x[2]**2,1)

nu=TrialFunction(V)

v=TestFunction(V)

a=dot(grad(nu), grad(v))*dx

L=rho*v*dx

nu=Function(V)

solve(a == L, nu, bc)

plot(nu)

#plot(mesh)
'''
