from formulas import *
from sympy.parsing.sympy_parser import parse_expr

def definedOrderMethods(params, formula):
    y0 = float(params[0])
    t0 = float(params[1])
    h = float(params[2])
    steps = float(params[3])
    f = parse_expr(params[4])

    i = 0
    yn = y0
    ypoints = [y0]
    while i < steps:
        yn1 = formula(f, yn, t0+(i*h), h)
        ypoints.append(yn1)
        yn = yn1
        i += 1
    
    return ypoints

def adamsBashforthMethod(params):
    order = int(params[len(params)-1])
    f = parse_expr(params[len(params)-2])
    steps = float(params[len(params)-3])-order+1
    h = float(params[len(params)-4])
    t0 = float(params[len(params)-5])
    prev_ypoints = params[0:order]

    j = 0
    coeffs = []
    while j < order:
        coeffs.append(calculateAdamsBashforthCoeffs(order, j))
        j += 1

    recurrence_relation = buildAdamsBashforthRecurrenceRelation(order, h, coeffs)

    i = 0
    while i < order:
        prev_ypoints[i] = float(prev_ypoints[i])
        i += 1
    
    i = 0
    ypoints = []
    ypoints.extend(prev_ypoints)
    while i < steps:
        ynk = adamsBashforthFormula(f, prev_ypoints, t0+(h*i), h, order, recurrence_relation)
        ypoints.append(ynk)
        prev_ypoints.append(ynk)
        prev_ypoints.pop(0)
        i += 1
    
    return ypoints

def adamsMoultonMethod(params):
    order = int(params[len(params)-1])-1
    f = parse_expr(params[len(params)-2])
    steps = float(params[len(params)-3])-order+1
    h = float(params[len(params)-4])
    t0 = float(params[len(params)-5])
    prev_ypoints = params[0:order]

    j = 0
    coeffs = []
    while j <= order:
        coeffs.append(calculateAdamsMoultonCoeffs(order, j))
        j += 1

    recurrence_relation = buildAdamsMoultonRecurrenceRelation(order, h, coeffs)

    i = 0
    while i < order:
        prev_ypoints[i] = float(prev_ypoints[i])
        i += 1
    
    i = 0
    ypoints = []
    ypoints.extend(prev_ypoints)
    while i < steps:
        ynk = adamsMoultonFormula(f, prev_ypoints, t0+(h*i), h, order, recurrence_relation)
        ypoints.append(ynk)
        prev_ypoints.append(ynk)
        prev_ypoints.pop(0)
        i += 1
    
    return ypoints

def backwardDifferentiationMethod(params):
    order = int(params[len(params)-1])-1
    f = parse_expr(params[len(params)-2])
    steps = float(params[len(params)-3])-order+1
    h = float(params[len(params)-4])
    t0 = float(params[len(params)-5])
    prev_ypoints = params[0:order]

    i = 0
    while i < order:
        prev_ypoints[i] = float(prev_ypoints[i])
        i += 1

    i = 0
    ypoints = []
    ypoints.extend(prev_ypoints)
    while i < steps:
        ynk = backwardDifferentiationFormula(f, prev_ypoints, t0+(h*i), h, order)
        ypoints.append(ynk)
        prev_ypoints.append(ynk)
        prev_ypoints.pop(0)
        i += 1
    
    return ypoints

def writePoints(params, ypoints, title):
    y0 = float(params[0])
    t0 = float(params[1])
    h = float(params[2])

    o = open('output.txt', 'a')

    o.write(title)
    o.write('y( {} ) = {}\n'.format(t0, y0))
    o.write('h = {}\n'.format(h))
    
    i = 0
    for yp in ypoints:
        o.write('{} {}\n'.format(i, yp))
        i += 1
    o.write('\n')

    o.close()

def euler(params):
    writePoints(params, definedOrderMethods(params, eulerFormula), 'Metodo de Euler\n')

def euler_inverso(params):
    writePoints(params, definedOrderMethods(params, backwardEulerFormula), 'Metodo de Euler Inverso\n')

def euler_aprimorado(params):
    writePoints(params, definedOrderMethods(params, improvedEulerFormula), 'Metodo de Euler Aprimorado\n')

def runge_kutta(params):
    writePoints(params, definedOrderMethods(params, rungeKuttaFormula), 'Metodo de Runge-Kutta\n')

def adam_bashforth(params):
    writePoints([params[0], params[len(params)-5], params[len(params)-4]], adamsBashforthMethod(params), 'Metodo de Adams-Bashforth\n')

def adam_bashforth_by_euler(params):
    euler_params = params[0:5]
    euler_params[3] = params[5]
    ypoints = definedOrderMethods(euler_params, eulerFormula)
    adam_params = ypoints + params[1:6]

    writePoints(params, adamsBashforthMethod(adam_params), 'Metodo de Adams-Bashforth por Euler\n')

def adam_bashforth_by_euler_inverso(params):
    backward_euler_params = params[0:5]
    backward_euler_params[3] = params[5]
    ypoints = definedOrderMethods(backward_euler_params, backwardEulerFormula)
    adam_params = ypoints + params[1:6]

    writePoints(params, adamsBashforthMethod(adam_params), 'Metodo de Adams-Bashforth por Euler Inverso\n')

def adam_bashforth_by_euler_aprimorado(params):
    improved_euler_params = params[0:5]
    improved_euler_params[3] = params[5]
    ypoints = definedOrderMethods(improved_euler_params, improvedEulerFormula)
    adam_params = ypoints + params[1:6]

    writePoints(params, adamsBashforthMethod(adam_params), 'Metodo de Adams-Bashforth por Euler Apromirado\n')

def adam_bashforth_by_runge_kutta(params):
    runge_kutta_params = params[0:5]
    runge_kutta_params[3] = params[5]
    ypoints = definedOrderMethods(runge_kutta_params, rungeKuttaFormula)
    adam_params = ypoints + params[1:6]

    writePoints(params, adamsBashforthMethod(adam_params), 'Metodo de Adams-Bashforth por Runge Kutta\n')

def adam_multon(params):
    writePoints([params[0], params[len(params)-5], params[len(params)-4]], adamsMoultonMethod(params), 'Metodo de Adams-Moulton\n')

def adam_multon_by_euler(params):
    euler_params = params[0:5]
    euler_params[3] = params[5]
    ypoints = definedOrderMethods(euler_params, eulerFormula)
    adam_params = ypoints + params[1:6]

    writePoints(params, adamsMoultonMethod(adam_params), 'Metodo de Adams-Moulton por Euler\n')

def adam_multon_by_euler_inverso(params):
    backward_euler_params = params[0:5]
    backward_euler_params[3] = params[5]
    ypoints = definedOrderMethods(backward_euler_params, backwardEulerFormula)
    adam_params = ypoints + params[1:6]

    writePoints(params, adamsMoultonMethod(adam_params), 'Metodo de Adams-Moulton por Euler Inverso\n')

def adam_multon_by_euler_aprimorado(params):
    improved_euler_params = params[0:5]
    improved_euler_params[3] = params[5]
    ypoints = definedOrderMethods(improved_euler_params, improvedEulerFormula)
    adam_params = ypoints + params[1:6]

    writePoints(params, adamsMoultonMethod(adam_params), 'Metodo de Adams-Moulton por Euler Aprimorado\n')
    
def adam_multon_by_runge_kutta(params):
    runge_kutta_params = params[0:5]
    runge_kutta_params[3] = params[5]
    ypoints = definedOrderMethods(runge_kutta_params, rungeKuttaFormula)
    adam_params = ypoints + params[1:6]

    writePoints(params, adamsMoultonMethod(adam_params), 'Metodo de Adams-Moulton por Runge Kutta\n')

def formula_inversa(params):
    writePoints([params[0], params[len(params)-5], params[len(params)-4]], backwardDifferentiationMethod(params), 'Metodo Formula Inversa de Diferenciacao\n')

def formula_inversa_by_euler(params):
    euler_params = params[0:5]
    euler_params[3] = params[5]
    ypoints = definedOrderMethods(euler_params, eulerFormula)
    bdf_params = ypoints + params[1:6]

    writePoints(params, backwardDifferentiationMethod(bdf_params), 'Metodo Formula Inversa de Diferenciacao por Euler\n')

def formula_inversa_by_euler_inverso(params):
    backward_euler_params = params[0:5]
    backward_euler_params[3] = params[5]
    ypoints = definedOrderMethods(backward_euler_params, backwardEulerFormula)
    bdf_params = ypoints + params[1:6]

    writePoints(params, backwardDifferentiationMethod(bdf_params), 'Metodo Formula Inversa de Diferenciacao por Euler Inverso\n')

def formula_inversa_by_euler_aprimorado(params):
    improved_euler_params = params[0:5]
    improved_euler_params[3] = params[5]
    ypoints = definedOrderMethods(improved_euler_params, improvedEulerFormula)
    bdf_params = ypoints + params[1:6]

    writePoints(params, backwardDifferentiationMethod(bdf_params), 'Metodo Formula Inversa de Diferenciacao por Euler Aprimorado\n')

def formula_inversa_by_runge_kutta(params):
    runge_kutta_params = params[0:5]
    runge_kutta_params[3] = params[5]
    ypoints = definedOrderMethods(runge_kutta_params, rungeKuttaFormula)
    bdf_params = ypoints + params[1:6]

    writePoints(params, backwardDifferentiationMethod(bdf_params), 'Metodo Formula Inversa de Diferenciacao por Runge Kutta\n')