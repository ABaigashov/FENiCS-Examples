
from fenics import *
from mshr import *
import numpy as np

R=50

R0=10

rho0=0.1

mesh=UnitSquareMesh(10, 10)

u_d=Expression('x[0]',degree=3)

V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition

def boundary(x, on_boundary):
	return on_boundary

bc = DirichletBC(V, u_d, boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

f = Expression('x[0]*x[0]<=100 ? 1 - (x[0]*x[0])/100 : 0', degree=3)

f=Constant(0)

f1=Expression('1/(x[0]+1)', degree=3)

f2=Expression('x[0]',degree=3)

b=f1*u.dx(1)*f1*v.dx(1)*dx*f2+u.dx(0)*v.dx(0)*dx*f2

L=f*v*dx


# Compute solution
u = Function(V)
solve(b == L, u, bc)

# Save to file and plot
File("classic.pvd") << u



'''
V = FunctionSpace(mesh, 'P', 2)

bc = DirichletBC(V, Constant(0), 'on_boundary')




nu=TrialFunction(V)

v=TestFunction(V)

a=dot(grad(nu), grad(v))*dx

L=rho*v*dx

nu=Function(V)

solve(a == L, nu, bc)

plot(nu)

#plot(mesh)
'''
