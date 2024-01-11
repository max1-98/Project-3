


from mpmath import *
import cmath
import math
mp.prec = 100
mp.dps = 100

# Exact Method
def ed1(x):

	return cos(pi*x)/(x)-sin(pi*x)/(pi*x**2)

def ed2(x):

	return -(pi*sin(pi*x)/x)+2*sin(pi*x)/(pi*x**3)-2*cos(pi*x)/(x**2)

def ed3(x):

	return ((3*pi**2*x**2-6)*sin(pi*x)+(6*pi-pi**3*x**3)*cos(pi*x))/(pi*x**4)

print("-------------------------")
print(ed1(1), "Exact First Derivative")
print(ed2(1), "Exact Second Derivative")
print(ed3(1), "Exact Third Derivative")
print("-------------------------")


# Complex Step Differentiation
h = 10**(-4)
x = complex(1,h)
y = complex(1,2*h)


def f(x):
	return sin(pi*x)/(pi*x)

print(f(x).imag/h, "CSD First Derivative")
print(-2*(f(x).real-f(1))/h**2, "CSD Second Derivative")
print(1/h**3*(2*f(x).imag-f(y).imag))
print("-------------------------")


# Forward Difference Method
h = 10**(-4)
print((f(1+h)-f(1))/h, "FD First Derivative")
print((f(1)+f(1+2*h)-2*f(1+h))/h**2, "FD Second Derivative")
print(3*(f(1+h)-f(1+2*h)+1/3*f(1+3*h)-1/3*f(1))/(h**3), "FD Third Derivative")
print("-------------------------")






