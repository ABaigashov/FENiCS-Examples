
from fenics import *
from mshr import *

R=50

R0=10

rho0=0.1

domain=Sphere(Point(0,0,0),R)

star=Sphere(Point(0,0,0),R0)

domain.set_subdomain(1,star)

mesh=generate_mesh(domain, 32)

plot(mesh)

markers=MeshFunction(’size_t’,mesh, 3, mesh.domains())

V = FunctionSpace(mesh, ’P’, 1)

bc = DirichletBC(V, Constant(0), ’on_boundary’)

dx = Measure(’dx’, domain=mesh, subdomain_data=markers)

rho=rho0*(1-x[0]**2-x[1]**2-x[2]**2)

nu=TrialFunction(V)

v=TestFunction(V)

a=dot(grad(nu), grad(v))*dx

L=rho*v*dx(1)

nu=Function(V)

solve(a == L, nu, bc)

plot(nu)

#plot(mesh)
