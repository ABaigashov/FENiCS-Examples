import numpy as np
from sympy import *
import math

f0=0.45

R0=5

x=Symbol('x')

N=1000

f=f0*(1-x*x/(R0*R0))


density = []
pressure = []
enthalpy = []

F = open('EOS.txt', 'r')

lines = F.read().split()

F.close()

for i in range(0, len(lines)-2, 3):
    enthalpy.append(float(lines[i]))
    density.append(float(lines[i+1]))
    pressure.append(float(lines[i+2]))

x0=2


def energy(X0):
    enthalpy0=f.subs(x,X0)
    for i in range(len(pressure)):
        for j in range(len(pressure)):
            if enthalpy0<enthalpy[j]:
                K1=(pressure[j]-pressure[j-1])/(math.exp(enthalpy[j])-math.exp(enthalpy[j-1]))
                K2=(density[j]-density[j-1])/(math.exp(enthalpy[j])-math.exp(enthalpy[j-1]))
                result1=pressure[j-1]+(math.exp(enthalpy0)-math.exp(enthalpy[j-1]))*K1
                result2=density[j-1]+(math.exp(enthalpy0)-math.exp(enthalpy[j-1]))*K2
                break
    return X0,result2,result1

G = open('Table.txt', 'w')
for i in range(N):
    X=i*R0/N
    RES=energy(X)
    a=[]
    a.append(RES[0])
    a.append(RES[1])
    a.append(RES[2])
    G.write(str(a[0]))
    G.write(' ')
    G.write(str(a[1]))
    G.write(' ')
    G.write(str(a[2]))
    G.write('\n')

G.close()
