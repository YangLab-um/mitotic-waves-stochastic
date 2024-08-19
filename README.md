A Julia script for stochastic simulation of cell cycle oscillations in 0D and 1D.
This supplements and cross-checks the simulation by the [main repository](https://github.com/DanielRuizReynes/mitotic_waves) of [our research](https://www.biorxiv.org/content/10.1101/2024.01.18.576267).

To run the Gillespie's SSA to generate realizations of Cdk1 dynamics time series, use `cell_cycle_2ode` from `src/gillespie.jl`.
The following example generates stochastic trajectories shown in Fig. S6.
```bash
cd mitotic-waves-stochastic
mkdir data
julia --project=. scripts/ssa.jl
```
The implemented cell cycle dynamics follows the model by [Yang el al. (2013)](https://pubmed.ncbi.nlm.nih.gov/23624406/).

Mitotic waves can be simulated using `cell_cycle_2ode` from `src/waves.jl`. Refer to the following script file.
```bash
julia --project=. scripts/langevin.jl
```
