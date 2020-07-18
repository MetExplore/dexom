# DEXOM
> Diversity-based enumeration of optimal context-specific metabolic networks

DEXOM is a Matlab library for the reconstruction and enumeration of diverse optimal context-specific metabolic networks. It requires COBRA Toolbox (included as submodule) and a MILP solver (CPLEX, Gurobi).

## Installation

Clone the project with the --recursive argument to clone also the COBRA Toolbox submodule.
```
git clone --recursive https://github.com/MetExplore/dexom.git
```

Once cloned, open the dexom folder in matlab and run the initialization script `dexomInit.m`. By default, it tries to use any COBRA Toolbox in the Matlab's path, and if COBRA is not detected, the embedded version is added to the path.

```
>> dexomInit;
DEXOM: (EX)traction of Optimal Metabolic-networks (v0.1.0)
This library was tested using Matlab 2015b, CPLEX v12.8 and COBRA Toolbox v3.0.6.

> Initializing DEXOM library for Matlab
 + Adding external dependencies...
	> Trying to initialize COBRA Toolbox ... Done (already loaded)
> IBM CPLEX selected as the default solver (v128)
> Testing DEXOM ... Done.
> DEXOM is ready to use.
```

During the initialization, DEXOM launches a few quick tests for network reconstruction and enumeration. After getting the previous message, the library is ready to use.

