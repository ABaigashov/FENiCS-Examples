from sympy import *

def approximator(f):
    """ Функция возвращает, аппроксимированную десятичными логарифмами,
        таблично-заданную функцию y = f(x), в терминах C++.

        На вход функции подается открытый для чтения файл формата .txt
        (f = open('table_func.txt', 'r'))
    """
    # Чтение и разбиение строк списка
    lines = f.read().split()
    f.close()

    # Создание списков для абсциссы и ординаты
    abscissa = []
    ordinate = []
    for i in range(0, len(lines)-1, 2):
        abscissa.append(float(lines[i]))
        ordinate.append(float(lines[i+1]))

    # Определение количества строк в списке
    N = len(abscissa)

    # Создание списка с неизвестными a_i параметрами для уравнения a*log(x + c)
    a_i = []
    for i in range(0, N, 1):
        a_i.append('a' + str(i))
        a_i[i] = symbols(a_i[i])

    # Создание системы аппроксимирующих функций с неизвестными параметрами
    # Где c = N - i
    x = symbols('x')
    f = 0
    for i in range(0, N, 1):
        a = a_i[i]*log(x + N - i)
        f = f + a

    # Определение системы уравнений с неизвестными параметрами
    eq = []
    for i in range(0, N, 1):
        eq.append(f.subs(x,abscissa[i]) - ordinate[i])

    # Определение затравочной вектор-функции для решения системы уравнений
    x0 = []
    for i in range(0, N, 1):
        x0.append(i)

    # Определение неизвестных параметров
    sol = solvers.nsolve((eq), (a_i), (x0))

    # Построение итоговой функции с определенными параметрами
    x = symbols('x[0]')
    func = 0
    for i in range(0, N, 1):
        a = sol[i]*log(x + N - i)
        func = func + a

    # Принтинг функции в язык С++
    f_code = printing.ccode(func)

    #print(f_code)

    return f_code

# f = open('table_func.txt', 'r')
# print(approximator(f))
