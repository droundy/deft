Dear editor,

...

Response to first referee

The first referee acknowledges that our paper is well written and that the field
is an important one, but does not believe that our work is a *'substantial
contribution*'' to the field.  The referee has three major criticisms of the
value of our work (aside from issues with the introduction section).

**We argue that our method is based on SAMC (as we state) and is a substantial
improvement on the method (convergence is of 1 - 2 orders of magnitude superior
as can be seen from Fig. 3 and Fig. 4)**

**In addition, we also provide a detailed comparison of 4 different methods on a
system that would potentially be of value to those most likely to implement our
method (such as the polymer community has implemented SAMC)**

**Furthermore, the guidlines that the referee is using come from Computational
Physics Papers in the Physical Review E (2013). This specifically states
'Rather, the requirement that "papers must contain new results in physics"
should be interpreted as a requirement of novelty with respect to the
method(s), the physics, or both.' As stated above it is only necessary that the
method is novel OR the physical problem is novel. The paper below accepted by
PRE confirms this thought.**

**Analysis of the convergence of the 1/t and Wang-Landau algorithms in the
calculation of multidimensional integrals** *(accepted by PRE)*
This paper compares the performance of two methods (WL and 1/t-WL) to calculate
different integrals (not a real physical system because the exact answer is
known).

Firstly, the referee criticizes our paper for testing the method on a single
physical system (the square well fluid).  The referee gives as an argument an
example of a paper which tests a method on three different physical systems.  We
give as counter examples *references that we cite* which each use a single
physical system as the test case.  We agree that our method should be tested on
additional different physical systems but dispute that this is necessary for the
first paper introducing this method.  

**Analysis of the convergence of the 1/t and Wang-Landau algorithms in the
calculation of multidimensional integrals** *(accepted by PRE)*
This paper compares the performance of two methods (WL and 1/t-WL) to calculate
different integrals (not a real physical system because the exact answer is
known).

**Avoiding boundary effects in Wang-Landau sampling** *(accepted by PRE)*
This paper with DP Landau as an author considers a single system (the 2D Ising
model). The novelty of the paper is providing a scheme to avoid boundary effects
in WL.

**Optimal modification factor and convergence of the Wang-Landau algorithm**
*(accepted by PRE)*
This paper examines the optimal modification factor for the WL algorithm and
tests on a single system (the 2D Ising model).

**Convergence of Stochastic Approximation Monte Carlo and modified Wang–Landau
algorithms: Tests for the Ising model**
*(accepted by Computer Physics Communications)*
This paper compares SAMC, WL, and 1/t-WL on a single system (the 2D Ising model).

**Stochastic approximation Monte Carlo and Wang–Landau Monte Carlo applied to a
continuum polymer model** *(accepted by Computer Physics Communications)*
SAMC and WL are compared for a single continuum model with a square well type
potential.

These papers help to highlight that the convergence properties of the various
Monte Carlo methods are of importance (while the system(s) are not necessarily).
In our case, we compare the convergence *four* different methods, one of which
(SAD) is an improvement on the existing algorithm (SAMC), on a single system
(square well fluid).

Secondly, the referee does not see much value in providing a temperature range
(rather than an energy range) as input to the method.  We believe that a
temperature range is *far* more physical and convenient, and this is the primary
value of our paper.  **We have modified the abstract to better highlight this
contribution.**

Thirdly, because we find that 1/t-WL converges at the same rate as our method,
the referee concludes that our method is not a substantial improvement. We note
first (as mentioned in the paper in *Section IV* that the 1/t-WL method is
solving a significantly easier problem, because it is told in advance the range
of energies that are of interest.  We agree, however, that the benefit of our
method lies specifically in the fact that a temperature range is given as an
input.  Achieving the same result as the SAD method with 1/t-WL would either
require multiple simulations over a wider range of energies, or would require
guessing in advance a wider range of energies than are needed.  In either case,
the convergence of 1/t-WL would be worse than our tests indicate.
*Jordan Says we could add* **We have modified the conclusion to better emphasize
this important distication**

**----------------------------------------------------------------------------**

1. The referee does not acknowledge that our method is based on SAMC (the
mathematical generalization of WL). In fact, the referee does not seem to
acknowledge *the existence of the SAMC method* which we find curious. We clearly
state in the paper in section I *'In this work, we have developed an improved
algorithm based on SAMC.'*  This is again reinforced in Secion III.

There are many authors that acknowledge SAMC as a mathematical generalization of
WL. We list a few here in our response.

**Multidimensional stochastic approximation Monte Carlo**
*(accepted by PRE)*
"the Wang-Landau MC [3] scheme and its mathematical generalization, stochastic
approximation MC (SAMC) [4–6].""

**Convergence estimation of flat-histogram algorithms based on simulation results**
*(accepted by Computer Physics Communications)*
"At the same time, Liang suggested a generalized simulation method:stochastic
approximation Monte Carlo (SAMC) [7,8] including the 1/t-WL as a particular case"

**Stochastic approximation Monte Carlo and Wang–Landau Monte Carlo applied to a
continuum polymer model**
*(accepted by Computer Physics Communications)*
"Liang et al.[30] showed that WLMC could be seen as a version of SAMC and using
the mathematical background of stochastic approximation methods they proved the
convergence of SAMC."

The referee further disputes the use of language *'we have introduced a new
algorithm'*. In response, we argue that this language is commonly used to
describe MC methods that have their foundation in other methods (SAD is a
derivative of SAMC). Our choice of language is further validated by other
researchers in the field that are published in PRE.

**Dynamically optimized Wang-Landau sampling with adaptive trial moves and
modification factors**
*(accepted by PRE)*
This work introduces an improved update scheme as the referee would call it.
The authors still use language such as *'To test our new algorithm...'*

Still we have changed the language to say *'we have introduced a new variant of
SAMC called SAD'"'* to better highlight the subtle distinction.

**----------------------------------------------------------------------------**

2. The referee suggests that we might include references to the originators of
the multicanonical methods.
**We have added a paragraph to Section I detailing the important impact of
multicanonical methods on influencing flat histogram methods. We have referenced
the original authors Berg and Neuhaus.**

**----------------------------------------------------------------------------**

3. The referee suggests that that it would be helpful to include references to
the Statistical *(not Stochastic)* Temperature Monte Carlo methods.
**We have added a paragraph to Section I detailing the important impact of STMC
and RESTMC and its improvement on the WL method and Replica Exchange
*respectively.**

**----------------------------------------------------------------------------**

4. The referee suggests that that it would be valuable to include references to
STMD and metadynamics.
**We have added a paragraph to Section I detailing the important impact of STMD
forging a connection between WL and metadynamics.**

**----------------------------------------------------------------------------**

5. The referee argues that our method should be tested on multiple other systems.
The referee gives for example a single paper accepted by *PRL* that compares STMC
to WL and introduces STMD.

**----------------------------------------------------------------------------**

Minor details:

6. The referee asks why we chose a maximum temperature of infinity.  We have
added a paragraph to Section III explaining the reasoning, which is that
infinite temperature allows the system to randomize more reliably, and costs us
essentially nothing.  We do note that for some systems or situations in which
the energy of maximum entropy is very high (or worse, infinite) we would need
to introduce a maximum temperature.

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

**----------------------------------------------------------------------------**

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
