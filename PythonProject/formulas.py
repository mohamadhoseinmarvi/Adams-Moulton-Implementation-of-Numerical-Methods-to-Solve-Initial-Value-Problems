from sympy import *
import math

def calculateAdamsBashforthCoeffs(order, j): 
    k = 0
    product = 1
    while k < order:
        if (k != j):
            product *= (Symbol('x')+order-k-1)
        k += 1
        
    return (((-1)**(order-j-1))/(math.factorial(j)*math.factorial(order-j-1)))*integrate(product, (Symbol('x'), 0, 1))

def calculateAdamsBashforthIntegral(order, coeffs):
    j = 0
    eq = 0
    for c in coeffs:
        eq += c*Symbol('fn+{}'.format(j))
        j += 1
    return eq

def buildAdamsBashforthRecurrenceRelation(order, h, coeffs):
    return Symbol('yn+{}'.format(order-1)) + h*calculateAdamsBashforthIntegral(order, coeffs)

def adamsBashforthFormula(f, ypoints, tn, h, order, recurrence_relation):
    j = 0
    for yp in ypoints:
        recurrence_relation = recurrence_relation.subs('fn+{}'.format(j), f.subs({'y': yp, 't': tn+(j*h)}))
        j += 1
        
    return recurrence_relation.subs('yn+{}'.format(order-1), ypoints[len(ypoints)-1])

def calculateAdamsMoultonCoeffs(order, j):
    k = 0
    product = 1
    while k <= order:
        if (k != j):
            product *= (Symbol('x')+order-k-1)
        k += 1

    return (((-1)**(order-j))/(math.factorial(j)*math.factorial(order-j)))*integrate(product, (Symbol('x'), 0, 1))

def calculateAdamsMoultonIntegral(order, coeffs):
    j = 0
    eq = 0
    for c in coeffs:
        eq += c*Symbol('fn+{}'.format(j))
        j += 1
    return eq

def buildAdamsMoultonRecurrenceRelation(order, h, coeffs):
    return Symbol('yn+{}'.format(order-1)) + h*calculateAdamsMoultonIntegral(order, coeffs)

def adamsMoultonFormula(f, ypoints, tn, h, order, recurrence_relation):
    j = 0
    for yp in ypoints:
        recurrence_relation = recurrence_relation.subs('fn+{}'.format(j), f.subs({'y': yp, 't': tn+(j*h)}))
        j += 1
    
    recurrence_relation = recurrence_relation.subs('fn+{}'.format(order), f.subs({'y': eulerFormula(f, ypoints[len(ypoints)-1], tn+((order-1)*h), h), 't': tn+(order*h)}))

    return recurrence_relation.subs('yn+{}'.format(order-1), ypoints[len(ypoints)-1])

def eulerFormula(f, yn, tn, h):
    fn = f.evalf(subs={'y': yn, 't': tn})

    return yn+(fn*h)

def backwardEulerFormula(f, yn, tn, h):
    fn1 = f.evalf(subs={'y': eulerFormula(f, yn, tn, h), 't': tn+h})

    return yn+(fn1*h)

def improvedEulerFormula(f, yn, tn, h):
    fn = f.evalf(subs={'y': yn, 't': tn})
    fn1 = f.evalf(subs={'y': eulerFormula(f, yn, tn, h), 't': tn+h})

    return yn + (fn+fn1)*(h/2)

def rungeKuttaFormula(f, yn, tn, h): # Considering 4th order
    kn1 = f.evalf(subs={'y': yn, 't': tn})
    kn2 = f.evalf(subs={'y': yn+((h*kn1)/2), 't': tn+(h/2)})
    kn3 = f.evalf(subs={'y': yn+((h*kn2)/2), 't': tn+(h/2)})
    kn4 = f.evalf(subs={'y': yn+(h*kn3), 't': tn+h})

    return yn + (kn1+(2*kn2)+(2*kn3)+kn4)*(h/6)

def backwardDifferentiationFormula(f, ypoints, tn, h, order):
    recurrence_relations = [
        Symbol('yn+0') + h*Symbol('fn+1'),
        (4/3)*Symbol('yn+1') - (1/3)*Symbol('yn+0') + (2/3)*h*Symbol('fn+2'),
        (18/11)*Symbol('yn+2') - (9/11)*Symbol('yn+1') + (2/11)*Symbol('yn+0') + (6/11)*h*Symbol('fn+3'),
        (48/25)*Symbol('yn+3') - (36/25)*Symbol('yn+2') + (16/25)*Symbol('yn+1') - (3/25)*Symbol('yn+0') + (12/25)*h*Symbol('fn+4'),
        (300/137)*Symbol('yn+4') - (300/137)*Symbol('yn+3') + (200/137)*Symbol('yn+2') - (75/137)*Symbol('yn+1') + (12/137)*Symbol('yn+0') + (60/137)*h*Symbol('fn+5'),
        (360/147)*Symbol('yn+5') - (450/147)*Symbol('yn+4') + (400/147)*Symbol('yn+3') - (225/147)*Symbol('yn+2') + (72/147)*Symbol('yn+1') - (10/147)*Symbol('yn+0') + (60/147)*h*Symbol('fn+6')
    ]

    recurrence_relation = recurrence_relations[order-1]

    j = 0
    for yp in ypoints:
        recurrence_relation = recurrence_relation.subs('yn+{}'.format(j), yp)
        j += 1
    
    recurrence_relation = recurrence_relation.subs('fn+{}'.format(order), f.subs({'y': eulerFormula(f, ypoints[len(ypoints)-1], tn+((order-1)*h), h), 't': tn+(order*h)}))

    return recurrence_relation
