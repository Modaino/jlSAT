# Boolean satisfyability with continuous time dynamical system
## Quick and easy implementation in julia
### SATCFC.jl 
This julia code implements the continuous time model for the Ising machine with SAT potential and chaotic feedback control
### noisySAT.jl 
This julia code implements a RODE model for the Ising machine with SAT potential, chaotic feedback control and quantum fluctuations (multiplicative noise model).
$$\frac{du}{dt}=f\left(u, p, t, \xi(t) \right)$$ where $\xi(t)$ is a random process.
### jlSAT.jl
This julia code implements the continuous time dynamical system (CTDS) model for solving SAT problem.
### SDESAT.jl
This julia code implements an SDE model for the Ising machine with SAT potential, chaotic feedback control and quantum fluctuations (additive noise model).
$$du=f(u,p,t)dt+g(u,p,t)d\xi$$
