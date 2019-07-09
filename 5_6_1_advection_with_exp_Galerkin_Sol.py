from fenics import *
import math

mesh = IntervalMesh(1000, 0.1, 1)

Q = FunctionSpace(mesh, 'P', 1)

def left(x, on_boundary):
    return on_boundary and x[0] == 0.1

import sympy as sym
x = sym.symbols('x[0]')
s = 1.0/x
s_code = sym.printing.ccode(s)
s = Expression(s_code, degree=1)
kappa = Constant(1.0/500.0)
f = Constant(0.0)

# ==============================================================================
# standard Galerkin solution :

leftBC = DirichletBC(Q, 1.0, left)
u = TrialFunction(Q)
v = TestFunction(Q)
us = Function(Q)
uf1 = Function(Q)

a = + kappa * u.dx(0) * v.dx(0) * dx + u * v * dx
L =  f * v * dx

solve(a == L, us, leftBC)

uf1.interpolate(us)

# ==============================================================================
# plotting :

from pylab import *

t = mesh.coordinates()[:, 0][::-1]
uf1 = uf1.vector().array()

fig = figure(figsize = (5, 3.5))
ax = fig.add_subplot(111)

ax.plot(t, uf1)

ax.set_xlim([0, 0.3])
tight_layout()
# savefig('reults/5_6_1_advection_with_exp_Galerkin_Sol.pdf')
savefig('reults/5_6_1_advection_with_exp_Galerkin_Sol.png')
