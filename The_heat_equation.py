from __future__ import print_function

#The important first line
from fenics import *

#Generating meshes and Defining the finite element function space
nx = ny = 8
mesh = UnitSquareMesh(nx, ny)
V = FunctionSpace(mesh, 'P', 1)


#Defining the boundary conditions
T = 2.0
num_steps = 10
dt = T/num_steps
alpha = 3
beta = 1.2 

u_D = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1]*x[0] + beta*t',
	degree=2, alpha=alpha, beta=beta,t=0)

def boundary(x, on_boundary):
	return on_boundary

bc = DirichletBC(V, u_D, boundary)
u_n = interpolate(u_D, V)


#Defining the trial and test functions
u = TrialFunction(V)
v = TestFunction(V)


#Defining the source term
f = Constant(beta - 2 - 2*alpha)


#Defining the variational problem
F = u*v*dx + dt*dot(grad(u), grad(v))*dx  - (u_n + dt*f)*v*dx
a, L = lhs(F), rhs(F)


#Forming and solving the linear system
u = Function(V)
t = 0
for n in range(num_steps):
	t += dt
	u_D.t = t
	solve(a == L, u, bc)
	plot(u)
	u_e = interpolate(u_D, V)
	u_n.assign(u)


# #--------------------------Ploting PVD------------------------
# #Plotting the solution using the plot command
# plot(u)
# plot(mesh)
#
# #Save solution to file in VTK format
# vtkfile = File('The heat equation/The heat equation.pvd')
# vtkfile << u



#--------------------------Ploting pylab------------------------
# get array componets and triangulation :
v = u.compute_vertex_values(mesh)
x = mesh.coordinates()[:,0]
y = mesh.coordinates()[:,1]
t = mesh.cells()

from pylab import *
from mpl_toolkits.axes_grid1 import make_axes_locatable

ax = plt.axes()

cm = get_cmap('viridis')
c = ax.tricontourf(x, y, t, v, 10, cmap = cm)
p = ax.triplot(x, y, t, '-', color='k', lw=0.2, alpha=0.4)

# Output in the file
for ext in ["png", "pdf"]:
    print("The_heat_equation.%s" % (ext,))
    plt.savefig("The_heat_equation.%s" % (ext,), bbox_inches="tight")
