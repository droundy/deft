Dear editor,

We appreciate the referees' careful reading of our paper, and feel
that their feedback has led to a significant improvement. We have
added new data and a section for the 31 atom Lennard-Jones cluster. We
compute the low-temperature specific heat for this test system and
show the convergence SAD vs the other methods quantitatively. A
quantitative comparison of convergence on the LJ31 system has not been
done in the literature to our knowledge.  We also added a note to the
Square-well fluid section highlighting that this system has not
previously been explored using broad histogram methods. We believe
that we have addressed all of the referees' concerns and per the first
referee's own words 'If that or similar data would be in the paper
there would be no doubt whatsoever that SAD is a "substantial
improvement of existing computational methods".'

...

Response to first referee
# -----------------------------------------------------------------------------#

In the previous review, the first referee acknowledged that  'The paper is very
well written. It is also true that the authors addressed all of my previous
concerns in a satisfactory manner.'

The first referee asked that we simulate (either 'faster and/or
reliably') low-temperature solid-solid transitions in specific
Lennard-Jones clusters. The referee gave LJ135 as an example, but
stated that if specific heat capacity curves were in the paper for
LJ135 or 'similar' then 'there would be no doubt whatsoever that SAD
is a "substantial improvement of existing computational methods".' We
chose to examine the LJ31 cluster for reasons we state in the paper,
specifically that this cluster has a clear solid-solid phase
transition at very low temperature that has been shown in the past to
be challenging to simulate using flat histogram methods.  We believe
that we have completely addressed the request by the first referee and
thank the referee for pointing our attention toward simulating
Lennard-Jones clusters.

We were surprised by how poorly the 1/t-WL algorithm performed on the
LJ31 system.  In particular, the convergence was very strongly
affected by the lower limit of the energy range chosen.  This example
thus strongly highlights the benefits SAD provides when examining the
low-temperature heat capacity, where the temperature of interest is
known, but the range of energies that affect that temperature are hard
to determine.

Response to second referee
# -----------------------------------------------------------------------------#
All of the second referee's comments have been addressed previously.
