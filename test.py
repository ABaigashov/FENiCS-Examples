import sympy as sym

x = sym.symbols('x[0]')

def f(a,b):
    return a*sym.log(x + b)

for i in range(0, 3, 1):
    for j in range(0, 3, 1):
        func = sym.sum(f(i, j))

print(func)
#f_code = sym.printing.ccode(f)

#print(f_code)
