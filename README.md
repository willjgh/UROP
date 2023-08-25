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

min k_{r}
s.t. \sum_{r=1}^{R} k_{r} * Q_{r} * P = 0
    k_{r} > 0

Given a sufficient number of equations (for the number of unknowns) the solution should be the exact k_{r} values used in simulation.

## Unkown P: bounds


