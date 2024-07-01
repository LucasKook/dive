# Distributional Instrumental Variables Estimation

The code in this repository reproduces the results in [1].

# Dependencies

This package depends on the `dare` package, which implements the loss function
and interfaces `deeptrafo` and `deepregression` for stochastic gradient descent
via `tensorflow` and `keras`. The `dare` package can be installed from 
[here](https://github.com/LucasKook/dare).

# Reproducibility

All results can be reproduced by running `make all` or executing the scripts in
`./inst/code/` manually following the order in the `Makefile`. The results can
also be reproduced in parts. For the 401k application, `make 401k`; for the
schooling application, `make schooling`; for the simulation results, `make
run-simulations vis-simulations`; Figure 3 can be reproduced with `make
loss-landscape`; and for all other figures use `make figures`.

# References

[1] Kook, L., & Pfister, N. (2024). Instrumental Variable Estimation of
Distributional Causal Effects. arXiv preprint arXiv:2406.19986.
[doi:10.48550/arXiv.2406.19986](https://doi.org/10.48550/arXiv.2406.19986).
