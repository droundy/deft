Dear editor,

We appreciate the referees' careful reading of our paper, and feel that their
feedback has led to a significant improvement.  We have made some significant
edits to the paper in response to their comments.  These edits were focused on
providing a broader background with regard to previous work, increasing clarity,
and general proofreading.  We believe that we have addressed all of the
referees' concerns with the exception of the unreasonable request by the first
referee that we present results on additional new systems.

...

Response to first referee

The first referee acknowledges that our paper is well written and that the field
is an important one, but does not believe that our work is a *'substantial
contribution'* to the field.  The referee has three major criticisms of the
value of our work (aside from issues with the introduction section).

Our method is based on SAMC (as we clearly state in the paper) and is a substantial
improvement on this method (convergence is of 1 - 2 orders of magnitude superior
as can be seen from Fig. 3 and Fig. 4).  In addition, we also provide a detailed
comparison with 4 different flat histogram methods on the square well system.

Firstly, the referee criticizes our paper for testing the method on a single
physical system (the square well fluid).  The referee gives as an argument an
example of a paper which tests a method on three different physical systems.  We
give as counter examples twelve of our citations [10-12, 15, 18, 19, 23, 35-39],
which are all papers (six of which were published in PRE) which examine flat
histogram methods using a single physical system. Many of these papers only
present the convergence properties of a single Monte Carlo method. In our case,
we compare the convergence of *four* different methods, one of which (SAD) is an
improvement on the existing algorithm (SAMC), on a single system (square well
fluid).  We also note that to our knowledge this work is the only instance of a
test of any flat histogram method to the square well fluid.  We agree that our
method should be tested on additional different physical systems but dispute
that this is necessary for the first paper introducing an improved method.

Secondly, the referee does not see much value in providing a temperature range
(rather than an energy range) as input to the method.  We believe that a
temperature range is *far* more physical and convenient, and this is one of the
primary contributions of our paper.  We have modified the abstract to better
highlight this contribution.

Thirdly, because we find that 1/t-WL converges at the same rate as our method,
the referee concludes that our method is not a substantial improvement. We note
first (as mentioned in the paper in *Section IV*) that the 1/t-WL method is
solving a significantly easier problem, because it is told in advance the range
of energies that are of interest.  We do agree that the benefit of our method
lies specifically in the fact that a temperature range is given as an input.
Achieving the same result as the SAD method with 1/t-WL would either require
multiple simulations over a wider range of energies, or would require guessing
in advance a wider range of energies than are needed.  In either case, the
convergence of 1/t-WL would be worse than our tests indicate. We have modified
the conclusion to better emphasize this important distinction.

**----------------------------------------------------------------------------**

1. The referee does not acknowledge that our method is based on SAMC (the
mathematical generalization of WL). In fact, the referee does not seem to
acknowledge *the existence of the SAMC method* which we find curious. The
referee states that "the original WL paper does not specify a scheme for
updating the update factor."  This is false.  Both of the 2001 WL papers specify
the same scheme that is generally used for the WL method (the update factor
drops by a factor of 2).  It is true that the papers do not *mandate* this
schedule.  However, neither our method nor SAMC can be considered an
implementation of the WL method, because the WL method *does* mandate that the
update factor is only decreased when the histogram has become sufficiently flat.

Our method clearly does *not* implement a new update scheme for the WL method.
We also do not simply implement a new update scheme for SAMC, since SAMC
specifies an update schedule which is predetermined.  As stated in Section I
*'In this work, we have developed an improved algorithm based on SAMC.'*  This
idea is again reinforced in Section III.

The referee further disputes the use of language *'we have introduced a new
algorithm'*. This language is commonly used to describe MC methods that have
their foundation in other methods (as SAD is a derivative of SAMC, and SAMC is
itself inspired by WL). However, we did revise the paper to clarify in the
conclusion that our algorithm is a variant of the SAMC method.

**----------------------------------------------------------------------------**

2. The referee suggests that we might include references to the originators of
the multicanonical methods.  We have added a paragraph to Section I detailing
the important impact of multicanonical methods on influencing flat histogram
methods. We have referenced the original authors Berg and Neuhaus.

**----------------------------------------------------------------------------**

3. The referee suggests that that it would be helpful to include references to
the Statistical *(not Stochastic)* Temperature Monte Carlo methods.  We have
added a paragraph to Section I detailing the important impact of STMC and RESTMC
and its improvement on the WL method and Replica Exchange respectively.

**----------------------------------------------------------------------------**

4. The referee suggests that that it would be valuable to include references to
STMD and metadynamics. We have added a paragraph to Section I detailing the
important impact of STMD forging a connection between WL and metadynamics.

**----------------------------------------------------------------------------**

5. The referee argues that our method should be tested on multiple other
systems. The referee gives for example a single paper accepted by *PRL* that
compares STMC to WL and introduces STMD.  The paper cited by the referee is
unusual in that it introduces *two* new methods (STMC and STMD), each of which
is valuable with a different type of system.  As we state above, we do not
believe this is a reasonable or ordinary expectation.  We cite six PRE papers
(as well as six other papers) that examine flat histogram methods using a single
system.

The referee mentions the possibility of mentioning ground states.  The goal of
our method is explicitly to *not* explore the entire energy range, which in most
cases will mean not either seeking or finding the ground state.

**----------------------------------------------------------------------------**

Minor details (according to the referee):

6. The referee asks why we chose a maximum temperature of infinity.  *We have
added a paragraph to Section III explaining the reasoning*, which is that
infinite temperature allows the system to randomize more reliably, and costs us
essentially nothing.  We do note that for some systems or situations in which
the energy of maximum entropy is very high (or worse, infinite) we would need
to introduce a maximum temperature.

7. The referee points out that on page 5 we say "SAD should perform
similarly..." and asks if there is any evidence? We have edited the paper to
make clear that we have not tested and confirmed the behavior of SAD in systems
with a continuum of energies.

8. The referee asks that we do more proofreading of the paper.  We have done so.
In the introduction, we improve the clarity of the paper by defining t as the
number of moves when first introducing 1/t-WL. As the first referee suggests we
edited phrases such as Monte Carlo and the update factor for consistency. We
corrected several typographical errors throughout the paper. References were
checked for completeness.

**----------------------------------------------------------------------------**

The second referee states that our paper "makes an important contribution to
improving algorithms for this type of system."  Furthermore, the second referee
recognizes the difficulty of setting the timescale parameter in SAMC, and
acknowledges that our algorithm addresses this challenge.  The referee did have
several specific questions, which we have numbered below.

**----------------------------------------------------------------------------**

1. The referee wrote that the description of the flat histogram methods in
section II was confusing, and should be rewritten to highlight the differences
between methods.  We have edited this section to better describe the
implementation of each method in more detail in section II, and to highlight the
major differences between the methods. The update factor is now shown in display
math for each method so that it can more easily be seen how each method differs.

**----------------------------------------------------------------------------**

2. The referee felt that the entropy argument in Section II.D was unclear, and
wished for the relationship between entropy and the update factor to be more
fully described, as well as the basis for Equation 7. We have added an
additional description of the relationship between entropy, ln w, and the update
factor immediately after equation 7 to make more clear how the total change in
entropy connects with the update factor schedule.

**----------------------------------------------------------------------------**

3. The referee asked for a better choice of line styles to highlight the SAD
method in our plots. We changed our line styles to make the comparison plots
easier to read, to be more effective when printed without color, and to
highlight the SAD method. We found that by thickening the SAD line style and
thinning the SAMC lines we could make it considerably easier to distinguish
between the different methods and to better highlight the SAD method in our
plots.

**----------------------------------------------------------------------------**

4. The referee finally suggested that it might be interesting to try out our
method with a Go protein model.  We agree that this could be interesting, and
will consider this suggestion for a future paper.
