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
