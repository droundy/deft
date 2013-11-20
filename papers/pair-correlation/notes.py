from __future__ import division
from sympy import *

sigma, delta, N, V, kT = symbols('sigma delta N V kT')
r1, r2, r3, r4, rN = symbols('r_1 r_2 r_3 r_4 r_N')
Theta = symbols('\\Theta', cls=Function)

allthetas = symbols(r'{\prod_{i<j}}\Theta')


def multint(expr, rs):
    for r in rs:
        expr = Integral(expr, r)
    return expr

Z = multint(allthetas, [r1, r2, r3, r4, rN])
Zsym = symbols('Z')


def n2(ra, rb):
    integrals = allthetas
    for r in [r1, r2, r3, r4, rN]:
        if r != ra and r != rb:
            integrals = Integral(integrals, r)
    return N**2 * integrals/Zsym

def n4(ra, rb, rc, rd):
    integrals = allthetas
    for r in [r1, r2, r3, r4, rN]:
        if r != ra and r != rb and r != rc and r != rd:
            integrals = Integral(integrals, r)
    return N**4 * integrals/Zsym

f = open('notes-content.tex', 'w')

f.write(r'''
To begin, let's define $Z$ as
\begin{equation}
%s
\end{equation}
''' % latex(Eq(Zsym, Z)))

n2sym = symbols('n_2',cls=Function)(r1,r2)
f.write(r'''
And now we define $n_2$ and $n_4$ (which we normally write as $n^{(2)}$
and $n^{(4)}$, but I will abbreviate in these notes):
\begin{equation}
%s
\end{equation}
''' % latex(Eq(n2sym, n2(r1, r2))))

n4sym = symbols('n_4',cls=Function)(r1,r2,r3,r4)
f.write(r'''\begin{equation}
%s
\end{equation}
''' % latex(Eq(n4sym, n4(r1, r2, r3, r4))))

Zplusdelta = Zsym + Rational(1,2)*multint(Zsym*n2(r1,r2)*Theta(r1,r2), [r1,r2]) + Rational(1,8)*multint(Zsym*n4(r1,r2,r3,r4)*Theta(r1,r2)*Theta(r3,r4), [r1,r2,r3,r4])

f.write(r'''
After some work that I won't show here, we find that $Z(\sigma+\delta)$ (referred
to here as $Z_+$) is given by:
\begin{equation}
%s
\end{equation}
''' % latex(Eq(symbols('Z_+'), Zplusdelta)))

Zplusdelta = Zplusdelta.subs(Zsym*n4(r1,r2,r3,r4), Zsym*n4sym)
Zplusdelta = Zplusdelta.subs(Zsym*n2(r1,r2), Zsym*n2sym)
f.write(r'''
Based on the definitions of $n_2$ and $n_4$, we have that
\begin{equation}
%s
\end{equation}
''' % latex(Eq(symbols('Z_+'), Zplusdelta)))


Zplusdelta = Zplusdelta.subs(n2(r1,r2), symbols('n_2',cls=Function)(r1,r2))

f.write(r'''
At this stage, we can substitute in the approximation that $n_4$ is a product of $n_2$s:
\begin{equation}
%s
\end{equation}
''' % latex(Eq(symbols('Z_+'), Zplusdelta)))

n = N/V
eta = pi/6*sigma**3*n

F = N*kT*eta*(4-3*eta)/(1-eta)**2

gsigma = diff(F, sigma)/(2*N*kT*n*sigma**2*pi)
gsigma = gsigma.subs(N, symbols('eta')*6/pi/sigma**3*V)
print latex(gsigma)
gsigma = simplify(gsigma)
print latex(gsigma)
gsigma = factor(gsigma)


f.write(r'''
\begin{equation}
%s
\end{equation}
''' % latex(Eq(symbols(r'g_\sigma'), gsigma)))

gsigmaprime = (diff(diff(F,sigma),sigma)/(N*kT*n*pi*sigma**2) - 2*gsigma/sigma)/2
gsigmaprime = gsigmaprime.subs(N, symbols('eta')*6/pi/sigma**3*V)
gsigmaprime = simplify(gsigmaprime)
gsigmaprime = factor(gsigmaprime)

f.write(r'''
\begin{equation}
  g_\sigma' = %s
\end{equation}
''' % latex(gsigmaprime))
