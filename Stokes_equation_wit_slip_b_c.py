from fenics import *
import matplotlib as mpl
import mshr

#-------------Create mesh and define function space-------------
domain = Rectangle(Point(0, 0), Point(1, 1))
mesh = generate_mesh(domain, 64)
sub_domains = MeshFunction('size_t', mesh, 0)
 
#-------------Taylor - Hood element-------------
V = VectorFunctionSpace(mesh, 'BDFM', 2)
Q = FunctionSpace(mesh, 'BDFM', 1)
W = V * Q

# variational problem :
u, p = TrialFunctions(W)
v, q = TestFunctions(W)

# no penetration boundary condition for velocity
# y = 0, y = 1 and around the dolphin :
u_n = Constant(0.0)

# inflow boundary condition for velocity at x = 1 :
u_O = Expression(("-sin(x[1]*pi)", "0.0"), degree=2)

# relavent measures :
ds = Measure("ds")[sub_domains]
dG_O = ds(0)
dG_r = ds(1)

# constants :
gamma = Constant(1e2)
h = CellSize(mesh)
n = FacetNormal(mesh)
I = Identity(2)
eta = Constant(1.0)
f = Constant((0.0,0.0))
beta = Constant(10.0)

def epsilon(u): return 0.5*(grad(u) + grad(u).T)
def sigma(u,p): return 2*eta*epsilon(u) - p*I

t = dot(sigma(u,p), n)
s = dot(sigma(v,q), n)

B_o = + inner(sigma(u,p),grad(v))*dx - div(u)*q*dx

B_g = - dot(n,t)*dot(v,n)*dG_O \
	- dot(u,n)*dot(s,n)*dG_O \
	+ gamma/h*dot(u,n)*dot(v,n)*dG_O \
	+ beta*dot(u,v)*dG_O \
	- inner(dot(sigma(u,p), n), v)*dG_r \
	- inner(dot(sigma(v,q), n), u)*dG_r \
	+ gamma/h*inner(v,u)*dg_r
F = + dot(f,v)*dx \
	+ gamma/h*U_n*dot(v,n)*dG_O \
	- iner(dot(sigma(v,q), n), u_O)*dG_r \
	+ gamma/h*inner(v,u_O)*dG_r

# solve variational problem
wh = Function(W)

solve(B_o + B_g == F, wh)

# get individual components with deep copy :
uO, u1 = uh.split(True)

#--------------------------Ploting------------------------
from pylab import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import colors, ticker

# calculate array componets :
vO = uO.compute_vertex_values(mesh)
v1 = u1.compute_vertex_values(mesh)
v = sqrt(vO**2 + v1**2 + 1e-16)
vO = vO/v
v1 = v1/v
x = mesh.coordinates()[:,0]
y = mesh.coordinates()[:,1]
t = mesh.cells()

ax = plt.axes()

v[v > 2.0] = 2.0
cm = get_camp('viridis')
ls = array([0.0,0.1,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.00001])
nm = colors.BoundaryNorm(ls, cm.N)
c = ax.tricontourf(x, y, t, v, camp=cm, norm=nm, levels=ls)
tp = ax.triplot(x, y, t, '-', color='k', lw=0.2, alpha=0.3)
q = ax.quiver(x, y, vO, v1, pivot='middle', color='k', scale=60, width=0.0015,
	headwidth=4.0, headlength=4.0, headaxislength=4.0)

#-----------Output in the file-------------------
for ext in ["png", "pdf"]:
    print("saving 3_3_Stokes_equation.%s" % (ext,))
    plt.savefig("3_3_Stokes_equation.%s" % (ext,), bbox_inches="tight")
