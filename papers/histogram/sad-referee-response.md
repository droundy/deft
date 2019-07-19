Dear editor,

...

Response to first referee

The first referee acknowledges that our paper is well written and that the field
is an important one, but does not believe that our work is a substantial
contribution to the field.  The referee has three major criticisms of the value
of our work (aside from issues with the introduction section).

Firstly, the referee criticizes our paper for testing the method on a single
physical system (the square well fluid).  The referee gives as an argument an
example of a paper which tests a method on three different physical systems.  We
give as counter examples references **look them up** (including references
**look them up** which the referee asked us to add) which each use a single
physical system as the test case.  We agree that our method should be tested on
additional different physical systems but dispute that this is necessary for the
first paper introducing this method.  **Decide whether to just throw in the
Ising...**

Secondly, the referee does not see much value in providing a temperature range
(rather than an energy range) as input to the method.  We believe that a
temperature range is *far* more physical and convenient, and this is the primary
value of our paper.  **We have modified the abstract to better highlight this
contribution.**

Thirdly, because we find that 1/t-WL converges at the same rate as our method,
the referee concludes that our method is not a substantial improvement. We note
first (as mentioned in the paper in **...**) that the 1/t-WL method is solving a
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

2. We have added a section detailing the important impact of multicanonical
methods on influencing methods like WL etc... referencing the original authors
to the introduction.

3.

4.

5.

Minor details:

6.

7.

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
