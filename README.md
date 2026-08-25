# Key rates for semi-DI QKD protocols

This repository contains Julia code for computing key rates of semi-DI QKD protocols.
The main ingredient is a generalization of the NPA hierarchy to the scenario where the operators on some subsystem are fixed. This is implemented in the directory `semi-NPA`.

The actual computation of the key rates is performed in [Pluto](https://plutojl.org/) notebooks. By convention, these files end in `_nb.jl`. For example, there are the notebooks: `min_ent_nb.jl` and `bff_nb.jl` which can be used to compute the entropy of one-sided BB84 protocols. To evaluate these notebooks, run the following terminal command at the root of this repository
```bash
julia --project=.
```
Then, in the Julia REPL, run
```julia
import Pluto; Pluto.run()
```
This opens a new browser tab with the Pluto interface. There, the relevant notebooks can be selected and evaluated.

## Documentation
The documentation for the NPA code can be found in `semi-NPA/docs/build/index.html`. Simply open this file in your web browser of choice.
