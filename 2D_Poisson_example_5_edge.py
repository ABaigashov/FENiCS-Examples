from fenics import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

#-------------Create DOLPHIN mesh and define function space----------
mesh = UnitSquareMesh (32, 32)
V = FunctionSpace (mesh, 'Lagrange', 1)
ff = FacetFunction ('size_t', mesh , 0)
V = FunctionSpace(mesh, 'CG', 2)

#------iterate through the facets and mark each if on a boundary :----------
#1 - West
# 2 - East
# 3 - North
# 4 - South
for f in facets(mesh):
	n = f.normal()
	tol = DOLFIN_EPS
	if n.x() <= -tol and n.y() < tol and f.exterior():
		ff[f] = 1
	elif n.x() >= tol and n.y() < tol and f.exterior():
		ff[f] = 2
	elif n.x() < tol and n.y() >= tol and f.exterior():
		ff[f] = 3
	elif n.x() < tol and n.y() <= -tol and f.exterior():
		ff[f] = 4

ds = Measure('ds')[ff]
dN = ds(3)
dS = ds(4)
dE = ds(2)
dW = ds(1)

#-------------Define boundary condition------------
u0 = Constant(1.0)
bcE = DirichletBC(V, u0, ff, 2)
bcW = DirichletBC(V, u0, ff, 1)
bc = [bcE, bcW]


#-------------Define variational problem------------
u = TrialFunction(V)
v = TestFunction(V)
f = Expression("100*(50*x[0]*x[0] + 60*x[1]*x[1])", degree=2)
g = Expression("x[1]*2", degree=2)
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*(dN + dS)

#---------------------Compute solution---------------------
u = Function(V)
solve(a == L, u, bc)

#--------------------------Ploting------------------------
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

#divider = make_axes_locattable(gca())
# cax = divider.append_axes('right', "5%", pad="3%")
# cbar = colorbar(c, cax=cax)

#-----------Output in the file-------------------
for ext in ["png"]:
    plt.savefig("results/3_1_Poisson_equation.%s" % (ext,), bbox_inches="tight")
