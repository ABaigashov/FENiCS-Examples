

from fenics import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
import math
from NumSrcF import numer_src_function

R0=0.0001

R1=30


tol=1E-14

#задание области, в которйй ищется решение

domain=Rectangle(Point(R0,0), Point(R1, np.pi/2))

#генерация сетки

mesh=generate_mesh(domain, 64)

#задание функционального пространства

SRC = FunctionSpace(mesh, 'CG', 1)  # so that dofs are only in mesh vertices

#количество элементов пространства (количество узлов, в которых определяется Функция
#источника)

n = SRC.dim()

print(n)

#размерность пространства (число измерений)

d = mesh.geometry().dim()

print(d)

#массив координат узлов
F_dof_coordinates = SRC.tabulate_dof_coordinates().reshape(n,d)
F_dof_coordinates.resize((n, d))

#массив x-координат узлов и y-ооординат узлов
F_dof_x = F_dof_coordinates[:, 0]
F_dof_y = F_dof_coordinates[:, 1]

#значение энтальпии в центре зведды (определяет значение плотности в центре)
g0=0.45

#приблизительный радиус звезды (для первого шага)
R0=5

x=sym.Symbol('x')
y=sym.Symbol('y')

"""функция, определяющая приблизительный профиль энтальпии (для первого шага)
Можно взять любую монотннную функцию, которая в нуле принимает значение g0, а при x=R0 равна нулю
"""

g=g0*(1-x*x/(R0*R0))

"""создание файла H.txt, в которой записываются построчно x и y-коорнннаты узлов
и значение энтальпии в этих узлах. Координаты узлов берутся из массивов F_dof_x, F_dof_y
"""

G = open('data/H.txt', 'w')
for j in range(n):
    X=F_dof_x[j]
    Y=F_dof_y[j]
    enthalpy0=g.subs(x,X)
    enthalpy0=enthalpy0.subs(y,Y)
    a=[]
    a.append(X)
    a.append(Y)
    a.append(enthalpy0)
    G.write(str(a[0]))
    G.write(' ')
    G.write(str(a[1]))
    G.write(' ')
    G.write(str(a[2]))
    G.write(' ')
    G.write('\n')

G.close()

"""задание имен файлов, которые содержат уравннние состояния (Name1) в формате
энтальпия - плотность энергии - давление и файла, содержащего значение энтальпии
в узлах координатной сетки (создан ранее)
"""

Name1='data/EOS.txt'
Name2='data/H.txt'

"""вызов функции которая дает массив значений плотности и давления в узлах сетки
"""
result=numer_src_function(Name1, Name2,n)

"""запись значений плотности и давления в узлах сетки в файл. Формат запиcи
x-координата узла - y-координата узла - плотность - давление
Для расчетов файл не требуется, запись просто для контроля и проверки
"""
G = open('data/Src.txt', 'w')
for j in range(n):
	a=[]
	X=F_dof_x[j]
	density=result[0][j]
	pressure=result[1][j]
	Y=F_dof_y[j]

	a.append(X)
	a.append(Y)
	a.append(density)
	a.append(pressure)
	G.write(str(a[0]))
	G.write(' ')
	G.write(str(a[1]))
	G.write(' ')
	G.write(str(a[2]))
	G.write(' ')
	G.write(str(a[3]))
	G.write('\n')

G.close()

"""Часть кода, связанная непосредственно с самой вариационной задачей
"""

#задание граничных условий для функции u_1 и u_2 на границе области
u_R=Expression(('0','0'),degree=2)

#задание граничных условий на нормальные производные на границе области
g=Constant(0.0)


P1 = FiniteElement('P', triangle, 2)
element = MixedElement([P1, P1])
V = FunctionSpace(mesh, element)

#пробные функции по числу неизвестных функций
v_1, v_2 = TestFunctions(V)

#неизвестные функции
u = Function(V)

u_1, u_2 = split(u)

#определение правой и левой границы. Для расчетов левая граница не требуется
def boundary_R(x, on_boundary):
	return on_boundary and near(x[0], R1, tol)

def boundary_L(x, on_boundary):
	return on_boundary and near(x[0], 0, tol)

#задание граничных условий Дирихле на правой граиице ("бесконечности")
bc_R = DirichletBC(V, u_R, boundary_R)

bcs=[bc_R]

'''
f = Expression('x[0]*x[0]<=0.25 ? 10*(0.25-x[0]*x[0]) : 0', degree=2)

f00=Expression('x[0]*x[0]<=0.25 ? 2*(0.25-x[0]*x[0]) : 0', degree=2)
'''

#задание плотности и давления в узлах сетки. Значения берутся из result выше
f = Function(SRC)

f00 = Function(SRC)

#для плотности считывается 1-ый столбец массива result, лля давления второй.
f.vector()[:] = result[0]
f00.vector()[:] = result[1]

"""f0 = множитель перед угловой компонентой градиента функции (1/r)
f1 - якобиан перехода в сферических координатах. Для задачи, где функции
зввисят от угла, необходимо будет добавить sin(x[1])
f2 - якобиан перехода в 4-мерных сферических координатах. Для задачи, где
функции зависят от угла, надо уточнить в справочнике!
f3 - якобиан перехода в полярных ооординатах (r). Для данного файла не требуется,
т.к. сферически симметричная задача сведена к системе двух уравнений, одно из
которых в 3-х мерной сферической с.к., другое - в 4-хмерной.
"""

f0=Expression('1/(x[0])', degree=2)

f1=Expression('x[0]*x[0]',degree=2)

f2=Expression('x[0]*x[0]*x[0]',degree=2)

f3=Expression('x[0]',degree=2)


"""Сама вариационная задача. Первое уравнение системы взято без изменения.
Второе для функции NA переписано дяя функции u_2=ln(NA)=nu+lnA=u_1+lnA, это достигается
путем деления обеих частей на NA и несложными преобразованиями. Третье уравнение
тогда не требуется.
Часть с лапласианом интегрируется по частям, остальные части остаются без изменений.
Каждое из уравнений умножается на тестовую функцию v_1 и v_2 соответственно.
Множитель A^2 выражен как exp(2*lnA)=exp(2*u_2-2*u_1)
"""

F=-f1*u_1.dx(0)*v_1.dx(0)*dx-f1*f0*f0*u_1.dx(1)*v_1.dx(1)*dx \
-f2*u_2.dx(0)*v_2.dx(0)*dx-f2*f0*f0*u_2.dx(1)*v_2.dx(1)*dx \
#скалярное произведение градиентов функций. Эта часть кода не требует доработки
#при переходе к задаче с вращеиием
-f1*v_1*u_1.dx(0)*u_2.dx(0)*dx-f1*f0*f0*v_1*u_1.dx(1)*u_2.dx(1)*dx \
-f2*v_2*u_2.dx(0)*u_2.dx(0)*dx-f2*f0*f0*v_2*u_2.dx(1)*u_2.dx(1)*dx \
#часть, связанная с нелинейными слагаемыми. Не требует доработки при перехоже к задаче
#с вращением
-f1*4*(np.pi)*(f+3*f00)*exp(2*u_2-1*1.3*u_1)*v_1*dx \
-f2*exp(2*u_2-2*u_1)*4*4*(np.pi)*f00*v_2*dx
#часть связанная с плотностью и давлением

#F=f2*u_1.dx(0)*v_1.dx(0)*dx+f2*f1*f1*u_1.dx(1)*v_1.dx(1)*dx \
#-f2*v_1*u_2.dx(0)*u_1.dx(0)*dx \
#+f3*u_2.dx(0)*v_2.dx(0)*dx+f3*f1*f1*u_2.dx(1)*v_2.dx(1)*dx \
#-f3*v_2*u_2.dx(0)*u_2.dx(0)*dx \
#-f2*f*exp(2*u_2-2*u_1)*v_1*dx-f3*exp(2*u_2-2*u_1)*f*v_2*dx

#запуск решения задачи
solve(F == 0, u, bcs)


#задание угловой координаты точки для построения радиального профиля функции
#при фиксированном угле
theta0=0

#списки
data_x=[]
data_y1=[]
data_y2=[]

#заполнение массивов для построения графиков
for i in range (1000):
    #радиальная коордиаата
	x1=0.0001+0.001*i*(R1)
	p=Point(x1,theta0)
    #значение функции u_1
	F1=u(p.x(), p.y())[0]
    #значение функции u_2
	F2=u(p.x(), p.y())[1]
    #запись в массивы
	data_x.append(x1)
	data_y1.append(F1)
	data_y2.append(F2)


#построение графиков

fig=plt.figure(figsize=(8,6), facecolor='pink', frameon=True)
plt.plot(data_x,data_y1)
plt.savefig('results/Profile1.png')
plt.show()

fig=plt.figure(figsize=(8,6), facecolor='pink', frameon=True)
plt.plot(data_x,data_y2)
plt.savefig('results/Profile2.png')
plt.show()


"""
Запуск следующих итераций решения задачи. Используя решение на 1-м шаге, необходимо
пересчитать энтальююю в узлах. Это деаается с использованием того, что
H(r)+nu(r)=const. Значение константы равно сумме g0 (изначально заданная энтальпия
в центре) и значению nu в центре. Таким образом, H=g0+nu(0)-nu(r),
"""

for k in range(8):
    #вычисление nu0 - значения nu=u_1 в центре
    p0=Point(0.0001,0)
    nu0=u(p0.x(), p0.y())[0]
    #вывод на экран необязателен, просто для контроля
    print(nu0)

    #запись пересчитннной энтальпии в узлах сетки в файл H1.txt.
    G = open('data/H1.txt', 'w')

    for i in range(n):
        a=[]
        p=Point(F_dof_x[i], F_dof_y[i])
        nu=u(p.x(), p.y())[0]
        #пересчет энтальпии H в данном узле
        H=g0+nu0-nu
        X=F_dof_x[i]
        Y=F_dof_y[i]
        a.append(H)
        G.write(str(X))
        G.write(' ')
        G.write(str(Y))
        G.write(' ')
        G.write(str(a[0]))
        G.write('\n')

    G.close()

    """повторение части кода, связанного с вариационной задачей. Заново объявляем
    функции v_1, v_2, u_1, u_2. Вычисляем плотность и давление в узлах сетки,
    используя файл H1 с пересчитанной энтальпией в узлах.
    Записываем нннкционал F. Решаем задачу. Строим графики.
    P.S. Граничнее условия не меняются, сетка не меняется.
    """
    v_1, v_2 = TestFunctions(V)

    u = Function(V)

    u_1, u_2 = split(u)

    Name1='data/EOS.txt'
    Name2='data/H1.txt'

    result=numer_src_function(Name1, Name2,n)

    f = Function(SRC)

    f00 = Function(SRC)

    f.vector()[:] = result[0]
    f00.vector()[:] = result[1]

    F=-f1*u_1.dx(0)*v_1.dx(0)*dx-f1*f0*f0*u_1.dx(1)*v_1.dx(1)*dx \
    +f1*v_1*u_1.dx(0)*u_2.dx(0)*dx \
    -f2*u_2.dx(0)*v_2.dx(0)*dx-f2*f0*f0*u_2.dx(1)*v_2.dx(1)*dx \
    +f2*v_2*u_2.dx(0)*u_2.dx(0)*dx \
    -f1*4*(np.pi)*(f+3*f00)*exp(2*u_2-2*u_1)*v_1*dx \
    -f2*exp(2*u_2-2*u_1)*4*4*(np.pi)*f00*v_2*dx

    solve(F == 0, u, bcs)

    theta0=0

    data_x=[]

    data_y1=[]

    data_y2=[]

    for i in range (1000):
    	x1=0.0001+0.001*i*(R1)
    	p=Point(x1,theta0)
    	F1=u(p.x(), p.y())[0]
    	F2=u(p.x(), p.y())[1]
    	data_x.append(x1)
    	data_y1.append(F1)
    	data_y2.append(F2)


    fig=plt.figure(figsize=(8,6), facecolor='pink', frameon=True)
    plt.plot(data_x,data_y1)
    plt.savefig('results/Profile1-1.png')
    plt.show()

    fig=plt.figure(figsize=(8,6), facecolor='pink', frameon=True)
    plt.plot(data_x,data_y2)
    plt.savefig('results/Profile2-1.png')
    plt.show()
