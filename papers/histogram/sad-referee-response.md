Dear editor,

...

Response to first referee

The first referee acknowledges that our paper is well written and that the field
is an important one, but does not believe that our work is a substantial
contribution to the field.  The referee has three major criticisms of the value
of our work (aside from issues with the introduction section).

***Jordan Says:*** Computational Physics Papers in Physical Review E
Appropriate papers for this section are those that report new advances in
computational physics. This includes papers that describe new developments of
computational methods or substantial improvements of existing computational
methods, provided they are applicable to significant physics problems.
In many cases such a paper will improve its impact by explicitly demonstrating
the usefulness of the described method, i.e., by applying it to some physics
problem (old or new). _However, it should be noted that such an application is
not a necessary condition for eligibility_. Rather, the requirement that "papers
must contain new results in physics" should be interpreted as a requirement of
novelty with respect to the method(s), the physics, or both.
*PRE only even loosely requires that new methods be tested on a single (old) or new system.
Simply arguing that our method is distinct from WL (the referees argument) 'it
is a new version of SAMC!' should be sufficient based on PRE guidelines.*

Firstly, the referee criticizes our paper for testing the method on a single
physical system (the square well fluid).  The referee gives as an argument an
example of a paper which tests a method on three different physical systems.  We
give as counter examples references **look them up**  which each use a single
physical system as the test case.  We agree that our method should be tested on
additional different physical systems but dispute that this is necessary for the
first paper introducing this method.  **Decide whether to just throw in the
Ising...**
It is quite common for papers detailing comparison of convergence for various
methods to only examine a single physical or sometimes even an unphysical (as in
the case of integrals) system. When examining the convergence and/or comparing
various computational methods it cannot therefore be expected that multiple
systems should be considered in the same work. The authors provide a list of
such works that they directly cite (many accepted by PRE) that prove their point:

Analysis of the convergence of the 1/t and Wang-Landau algorithms in the calculation
of multidimensional integrals *(accepted by PRE)*
This paper compares WL and 1/t-WL on several different integrals (not even a
real physical system because the exact answer is known)

Avoiding boundary effects in Wang-Landau sampling *(accepted by PRE)*
This paper with DP Landau as an author only considers the 2D ising model with
the novelty of the paper to avoid boundary effects in WL which is succinctly
stated in the abstract: 'A simple modification of the ‘‘Wang-Landau sampling’’
algorithm...'

Optimal modification factor and convergence of the Wang-Landau algorithm *(accepted by PRE)*
This paper examines the optimal modification factor and tests on a single
system: the 2D ising model.

Convergence of Stochastic Approximation Monte Carlo and modified
Wang–Landau algorithms: Tests for the Ising model *(accepted by Computer Physics Communications)*
This paper compares SAMC, WL, and 1/t-WL on a single system: the 2D ising model

Stochastic approximation Monte Carlo and Wang–Landau Monte Carlo applied to a continuum polymer model *(accepted by Computer Physics Communications)*
SAMC and WL are compared for a single continuum model with a square well type potential.

Secondly, the referee does not see much value in providing a temperature range
(rather than an energy range) as input to the method.  We believe that a
temperature range is *far* more physical and convenient, and this is the primary
value of our paper.  **We have modified the abstract to better highlight this
contribution.**

Thirdly, because we find that 1/t-WL converges at the same rate as our method,
the referee concludes that our method is not a substantial improvement. We note
first (as mentioned in the paper in **...** The results section we could highlight again in the conclusion?) that the 1/t-WL method is solving a
significantly easier problem, because it is told in advance the range of
energies that are of interest.  We agree, however, that the benefit of our
method lies specifically in the fact that a temperature range is given as an
input.  Achieving the same result as the SAD method with 1/t-WL would either
require multiple simulations over a wider range of energies, or would require
guessing in advance a wider range of energies than are needed.  In either case,
the convergence would be worse than our tests indicate.

1. **Existence of the SAMC method.** We state in the paper "In this work, we
have developed an improved algorithm based on SAMC."  This is again reinforced
in Secion III.

Multidimensional stochastic approximation Monte Carlo
*(accepted by PRE)*
"the Wang-Landau MC [3] scheme and its mathematical generalization, stochastic approximation
MC (SAMC) [4–6].""

Convergence estimation of flat-histogram algorithms based on simulation results
*(accepted by Computer Physics Communications)*
"At the same time, Liang suggested a generalized simulation method:stochastic
approximation Monte Carlo (SAMC) [7,8] including the 1/t-WL as a particular case"

Stochastic approximation Monte Carlo and Wang–Landau Monte Carlo applied to a continuum polymer model
*(accepted by Computer Physics Communications)*
"Liang et al.[30] showed that WLMC could be seen as a version of SAMC and using
the mathematical background of stochastic approximation methods they proved the
convergence of SAMC."

2. We have added a section detailing the important impact of multicanonical
methods on influencing methods like WL etc... referencing the original authors
to the introduction.

3. As per the referees request, we have briefly reviewed methods the Statistical
temperature Monte Carlo method.

4.  As per the referees request, we have briefly reviewed methods used in the
metadynamics communities in the introduction.

5. more test cases

Minor details:

6. Out of curiosity: The authors say a temperature range T_min < T <
infinity should be chosen. Why not T_min < T < T_max?

7. On page 5 the authors say "SAD should perform similarly ...". Is
there any evidence?

8. *Proofreading* **Note all changes made under proof reading, as well as what was
suggested.**
In the introduction, we improve the clarity of the paper by defining t as the number
of moves when first introducing 1/t-WL.
As the first referee suggests we edited phrases such as Monte Carlo and the
update factor for consistency.
We corrected several typographical errors throughout the paper.
References were checked for completeness.



The second referee was far more positive on the value of our work, but did
suggest several places where the paper could be made more clear.

1. The referee wrote that the description of the flat histogram methods in
section II was confusing, and should be rewritten to highlight the differences
between methods.  We have added ...

2. The referee shared that the entropy argument in Section II.D is unclear, and
wished for the relationship between entropy and the update factor to be more
fully described, as well as the basis for Equation 7.

3. The referee asked for a better choice of line styles to highlight the SAD
method in our plots.

4. The referee finally suggested that it might be interesting to try out our
method with a Go protein model.  We agree that this would be interesting, and
will consider this suggestion for a future paper.
