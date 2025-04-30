# stochasticmodellingcoursework
Stochastic Modelling Coursework 2024

This project models the dynamics of extra-chromosomal circular DNA, or ecDNA, in a population of dividing cells. In this model, we denote by Nk the number of cells that contain k copies of ecDNA, where k is a non-negative integer.

The model is governed by two key differential equations:

The first equation describes the rate of change in the number of cells that have zero copies of ecDNA. This rate increases in two ways: first, due to natural cell division of cells that already have no ecDNA, which we assume happens at a constant rate of one; and second, due to cells that previously had ecDNA but lost all copies during cell division. The latter process is modelled using a probability factor that accounts for the chance that all ecDNA copies go to just one daughter cell during division. Since each dividing cell produces two daughter cells at a certain division rate, the equation includes a factor of two times the division rate.

The second equation describes the rate of change in the number of cells with a non-zero number of ecDNA copies. This rate decreases due to the division of cells with exactly k copies, as these cells are being replaced. It also increases due to other cells that divide and generate daughter cells with exactly k copies of ecDNA. The increase is calculated by summing over all parent cells that have enough ecDNA to produce such a daughter cell, specifically, cells with at least k divided by two copies. For each possible number of ecDNA copies in a parent cell, the equation considers both the number of ways that k copies can be distributed among the total available copies, and the probability of that particular distribution occurring.

Together, these equations capture how ecDNA is inherited or lost through stochastic partitioning during cell division, providing insight into how ecDNA dynamics evolve over time in a growing cell population.

