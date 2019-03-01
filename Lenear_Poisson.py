from __future__ import print_function
from fenics import *

# Create mesh and define function space
mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition
u_D = Expression('1 + x[0]*x[0] + x[1]*x[1]', degree=2)

def boundary(x, on_boundary):
	return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u), grad(v))*dx
L = f*v*dx

# Compute solution
u = Function(V)
solve(a == L, u, bc)


# #--------------------------Ploting PVD------------------------
# # Plot solution and mesh
# plot(u)
# plot(mesh)
#
# # Save solution to file in VTK format
# vtkfile = File('Lenear_Poisson/Lenear_Poisson.pvd')
# vtkfile << u
#
# # Compute maximum error at vertices
# vertex_values_u_D = u_D.compute_vertex_values(mesh)
# vertex_values_u = u.compute_vertex_values(mesh)



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
    print("Lenear_Poisson.%s" % (ext,))
    plt.savefig("Lenear_Poisson.%s" % (ext,), bbox_inches="tight")
