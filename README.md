# UROP
Code written for a UROP project on Markov chain inference in cell biology

# Reaction Test
Initial investigation into simulating trajectories of reaction chains using the Gillespie algorithm.

## Gillespie

Given an average of \lambda independent arrivals per unit time, the number of arrivals per unit time is X ~ Poisson(\lambda), and the number of arrivals in the next t units of time is X_{t} ~ Poisson(\lambda t). This gives a "waiting time" until the next arrival of T ~ Exp(\lambda)

Given the chain is in state i, it will stay there for an Exp(-qii) distributed "holding time", then transitioning to a state j with probability - qij / qii. To simulate a trajectory we simply need to simulate these holding times and transitions.

Steps:
- initial state i
- simulate Exp(-qii) time
- choose state j to transition to with prob - qij / qii
- repeat for state j, until time limit reached

## Stationary distributions
Simulating a sample path/trajectory for a sufficient length of time (depending on transition rates and context), we can asssume that the chain has entered its stationary regime. The final state (at stopping time) is a sample from the stationary distribution of the chain. So a simple stationary estimate is given by sampling many trajectories and calculating the proportion of each state among the final states.

However, this is just a point estimate with no error guarantees. A better estimate would give confidence intervals for the stationary distribution value of each state:

### Bootstrapping
- Simulate sample path up to sufficient time
- record final state: sample from stationary distribution
- repeat for n samples
- Bootstrap: sample with replacement to obtain N bootstrap samples, each of size n
- calculate point estimate of stationary distribution for each bootstrap sample (proportions of selected states within the sample)
- Confidence Intervals: for each state have N estimates of stationary probaility, use e.g. 97.5% and 2.5% quantiles for 95% CI for true prob.

Plots can be produced to show histograms of stationary estimates along with CI bounds, for each state. 

Introductory analysis of variance using width of CI is shown, but expanded later.

# Linear Prog
Initial investigation into using packages to solve LP optimization problems, and formulating Qp = 0 systems to solve for parameters.

## LP
The true statinonary distribution p satisfies Qp = 0, where Q is the rate matrix (or generator G) of the chain/reaction system containing the reaction rate constants k_{r}: parameters that we wish to infer. This can be decomposed into a sum of constant Q_{r} matrices:
Qp = \sum_{r=1}^{R} k_{r} * Q_{r} * P = 0
#### Outer Approximation
The true reaction rate constants k_{r} must satisfy this equation (define as membership of a set S), which can then be used as a constraint for LP to find upper and lower bounds via min and max of k_{r}. The solution bounds provide an "outer approximation": a larger, more conservative set, containing the true set of solutions, which can then be refined, hopefully approcahing the true set of k_{r} values.
### Infinite equations
The state space of species in reaction (e.g. number of molecules) is generally infinite (0,1,2,...), giving an infinite system of equations (infinite dim Q). However, we cannot make use of all of them, nor do we need to, using *state space truncations*.
In this, we choose a subset of states x (or equations/rows of Q) to use, only including equations involving those P(x). Initially we use the first few x's: 0,1,2,... but we will later investigate what the best truncations to use are, and how this depends on context.

## Known P
Assuming the exact stationary distribution P is known, the main LP equation reduces to a set of linear equality constraints, with a linear objective, for the k_{r} variables:

min/max k_{r}
s.t. \sum_{r=1}^{R} k_{r} * Q_{r} * P = 0
    k_{r} > 0

Given a sufficient number of equations (for the number of unknowns) the solution should be the exact k_{r} values used in simulation.

## Unkown P: bounds
If we instead have bounds, e.g. bootstrap CI's, on the true values p(x) then we must include them as variables in the LP: P_{l} < P < P_{u}. However, this introduces non-linear terms into the Qp = 0 constraint: namely k_{r} * P, cross terms of k_{r} * p(x). To linearise we introduce z variables: z_{r} = k_{r} * P, and write the constraints in terms of z's and k's:

min/max k_{r}
s.t. \sum_{r=1}^{R} z_{r} * Q_{r} = 0
    k_{r} * P_{l} < z_{r} < k_{r} * P_{u}
    k_{r} > 0

(Note: this is possible because only cross terms are present, there are no quadratic terms e.g. p(x)^2)

## Packages
PuLP was inconvenient, but CVXPY has very nice vector and matrix variable support, and makes it easy to solve both the min and max problem and report the outcome (optimal/unbounded/infeasible). The code shows birth-death (M/M/1 queue: death rate = 1) and schlogl examples.

# Birth Death Reaction
Investigate a simple birth-death process of one species, look at relationship of variance in estimates with true value, and the best state space truncations to use.

## Reaction

\empty -{k1}-> M

M -{k2}-> \empty

Molecules of M are produced with a rate of k1 and decay a rate of k2 (per molecule!). The stationary distribution is Poisson(k1 / k2) (e.g. p(x) = P(X = k) where X ~ Poi(k1/k2)). 

Code to simulate sample paths of the reaction, and to find bootstrap confidence intervals for stationary distribution values is provided.

## Confidence Intervals: width v true value

The width of each bootstrap CI is a measure of the variance/sampling error of the estimate of the stationary distribution p(x), but how does this relate to the true value of p(x)?

Plotting CI width against p(x) for a selection of states x and a range of parameter values (k1, k2) shows that width \prop p(x) * (1 - p(x)). However, the raw width is a measure of abosolute error: when sampling even 1000 final states for a state of prob 0.5 we would expect a deviation of around 450-550, so the CI width could vary by ~0.1, but for a state of prob 0.005 even with a large variation in the sample e.g. 0 to 10 the CI width would only vary by ~0.01. So we also consider the relative error: CI width / true which decreases as p(x) increases. This is perhaps more intuitive; as the true prob increases to 1, we observe the state more often and so have more information, whereas as the prob decreases to 0 we observe it less and less, giving less information with which to estimate (in the extreme case of p(x) = 0, we will never observe the state, so we estimate 0 but with no degree of uncertainty, and if the true value was actually e.g. 0.0001 we would also be unsure). However, we see the relative error is \prop 1 / p, with a rapid decrease from 0 to ~0.1, levelling off as sufficient events are observed to estimate well.


## State Space Truncations: What are the best equations to use?
The main constraint of the LP is the equation \sum_{r=1}^{R} Q_{r} * z_{r} = 0, where Q_{r} are constant matrices of bi/tri/etc-diagonal structure. The nth equation/row of this system involes p(n - 1), p(n) and p(n + 1) (due to upper and lower bidiagonal Q_{r}'s). Since the width of CI's for the estimates of these values varies, we can ask which equations, and so which states will give the tightest or most accurate solution bounds. Of course in practice we can use as many equations as possible, but some states may have confidence intervals that do not include the true value (as seen in the previous section with low prob states having greater variability) leading to infeasible LPs, so finding the best states to use can help avoid those that may cause problems.

Code is provided to solve birth death LPs given upper and lower confidence bounds, as well as a list of equations to use (note: the nth equation involves p(n-1),p(n),p(n+1) as seen, so the bounds must be provided up to p(n+1)).

### Best Single Equation
Solving the LP using only a single equation/row of Qp = 0: simulate CI bounds, find solution interval using each equation alone (up to a given value) and plot: plot width vs equation, and also visual intervals.

See that the first equations, espcially 1st, seem to give the best bounds, even for distributions centered around higher states than appear in them. Also see that the stationary distribution (Poi(k1/k2)) can be bimodal for integer parameters, leading to no upper bounds on equations involving modal states, even if they should contain the most information (p(x)'s equal so cancel out???).

To find the best equations repeat this many times, recording the equation of minimum width, then plot the distribution of best equations.

The distributions produced mostly peak at 0: likely due to the fewer terms in the 1st equation, as well as the non-negative number of molecules leading to higher p(0) values, reducing the CI width as seen before. However, even with stationary distributions centered around e.g. 1, the best equation distribution produces small peaks at e.g. equations 10 - 12, which is unexpected as these equations are not expected to contain information.

### Multiple equations
Given simulated data, solve the LP using a single equation and find the solution interval then repeat, adding an extra row/equation each time. Code is provided to do this starting with the 1st equation adding until the nth, and the nth equation adding until the first.

In general there is a sharp decrease in solution width when adding the first few equations, then a shallow decrease until the equation containing the mode states of the distribution, where the solution width stagnates and does not improve with additional information. The same is broadly true of adding equations in descending order, but the 1st equation still seems to provide extra information, again perhaps due to fewer terms or a smaller CI.

# Gene Expression
Code exploring simulation and solution of a model for gene expression.

# Non-Linear Programming
So far we have considered only linear programs, formed from the reaction systems by using z_{r} = k_{r} * P variables to linearise the bilinear constraints, but we can also directly use these bilinear constraints to form a bilinear program. However, these constraints are non-convex (Q_{r} matrices are not symmetric, let alone PSD etc) and so we cannot use CVXPY which only deals with convex problems and instead work with scipy.optimize. Since we have lost convexity there is now no guarantee of global optimality, and we must be careful to avoid becoming stuck in local optima.

This approach works well on the birth-death reaction, giving mostly tight and accurate bounds on the birth rate constant, although sometimes the bounds are innacurate and do not contain the true value.

However, when attempting the gene expression model the results are poor with most solution bounds very far from the true values, perhaps due to a more complex landscape and the prescence of many local minima and maxima.




