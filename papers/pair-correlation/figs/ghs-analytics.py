from __future__ import division
import sympy
from sympy import pi, exp, symbols

eta_s = symbols(r'\eta')
eta = eta_s
K0 = symbols(r'\kappa_0')
K1 = symbols(r'\kappa_1')
K2 = symbols(r'\kappa_2')

r = symbols(r'r')
sigma_s = symbols(r'\sigma')
sigma = sigma_s
R_s = symbols(r'R')
R = sigma/2

z_s = symbols('z')
z = (r - sigma)/R

hsigma = (1 - eta/2)/(1 - eta)**3 - 1
hsigma_equation = sympy.Eq(symbols(r'h_\sigma'), hsigma)
#all_etas = sympy.solve(hsigma_equation, eta)
#eta = all_etas[0]
hsigma_s = symbols(r'h_\sigma')
density_s = symbols(r'n')
density = 3/(4*pi)*eta
density = density_s

rhs = (1-eta)**4/(1 + 4*eta + 4*eta**2 - 4*eta**3 + eta**4)

int_h0 = 4*pi*hsigma*(2 + sigma*K0*(2 + sigma*K0))/K0**3 - 4*pi*sigma**3/3

int_h1_over_a1 = 4*pi*(6 + sigma*K1*(4 + sigma*K1))/K1**4

int_h2_over_a3 = 8*pi*(12 + sigma*K2*(6 + sigma*K2))/K2**5

A = (rhs-1)/density - int_h0
B = int_h2_over_a3
C = int_h1_over_a1

a1_s = symbols(r'a_1')
a1 = hsigma*(K0 - 1 - hsigma) # sets slope at sigma

a3_s = symbols(r'a_3')
a3 = (A - C*a1)/B # sets integral

print r"""
\documentclass{article}

\usepackage{breqn}

\begin{document}"""

hsigma = hsigma_s
z = z_s
a1 = a1_s
a3 = a3_s
ghs = 1 + hsigma*exp(-K0*z) + a1*z*exp(-K1*z) + a3*z**2*exp(-K2*z)

print r"""
\begin{dmath}
  A = %s
\end{dmath}
""" % sympy.latex(A.subs(z, symbols('z')))

print r"""
\begin{dmath}
  B = %s
\end{dmath}
""" % sympy.latex(B.subs(z, symbols('z')))

print r"""
\begin{dmath}
  a_3 = %s
\end{dmath}
""" % sympy.latex(a3.subs(z, symbols('z')))

print r"""
\begin{dmath}
  g_{HS}(r) = %s
\end{dmath}
""" % sympy.latex(ghs) # sympy.latex(sympy.simplify(ghs))


print r"""
\end{document}
"""
