from sympy import *

def approximator(f):
    """ Функция возвращает, аппроксимированную десятичными логарифмами,
        таблично-заданную функцию, в терминах C++
    """
    # Чтение списка
    lines = f.read().split()
    f.close()

    # Создание списков для абсциссы и ординаты
    abscissa = []
    ordinate = []
    for i in range(0, len(lines)-1, 2):
        abscissa.append(float(lines[i]))
        ordinate.append(float(lines[i+1]))

    # Количество строк в списке
    N = len(abscissa)

    # Создание списка с параметрами для уравнения a*log(x + c)
    k_a = []
    for i in range(0, N, 1):
        k_a.append('a' + str(i))
        k_a[i] = symbols(k_a[i])

    # Создание системы аппроксимирующих функций с неизвестными параметрами
    x = symbols('x')
    f = 0

    for i in range(0, N, 1):
        k = k_a[i]*log(x + N - i)
        f = f + k

    # Определение системы уравнений с неизвестными параметрами
    eq = []
    for i in range(0, N, 1):
        eq.append(f.subs(x,abscissa[i]) - ordinate[i])

    # Определение затравочной вектор-функции для решения системы уравнений
    x0 = []
    for i in range(0, N, 1):
        x0.append(i)

    # Определение неизвестных коэффициентов
    sol = solvers.nsolve((eq), (k_a), (x0))

    # Построение итоговой функции с определенными коэффициентами
    x = symbols('x[0]')
    func = 0
    for i in range(0, N, 1):
        k = sol[i]*log(x + N - i)
        func = func + k

    # Принтинг функции в язык С++
    f_code = printing.ccode(func)

    print(f_code)

    return f_code

#approximator()
