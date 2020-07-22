# DEXOM
> Diversity-based enumeration of optimal context-specific metabolic networks

DEXOM is a Matlab library for the reconstruction and enumeration of diverse optimal context-specific metabolic networks. It requires COBRA Toolbox (included as submodule) and a MILP solver (CPLEX, Gurobi).

<p align="center"><img src="https://github.com/MetExplore/dexom/raw/master/assets/overview.png" width="500"></p>

## Installation

In order to use DEXOM, you need a compatible version of Matlab (recommended version >= 2015), CPLEX (recommended version >= 12.8), and COBRA Toolbox (recommended v3.0.6). This repository includes an embedded version of COBRA Toolbox v3.0.6 with minor changes to facilitate the reproducibility of the experiments and avoid automatic updates that can introduce incompatibilities. In order to clone DEXOM with the embedded COBRA v3.0.6, use the following git command:

```
git clone --recurse-submodules="modules/cobratoolbox" https://github.com/MetExplore/dexom.git
```

Once cloned, open the dexom folder in matlab and run the initialization script `dexomInit.m`. By default, it tries to use any COBRA Toolbox in the Matlab's path, and if COBRA is not detected, the embedded version is added to the path. In order to prevent incompatibilities with your current COBRA version, it is recommended to remove all COBRA folders from your Matlab's path before initializing DEXOM, in order to force DEXOM to use the embedded COBRA Toolbox.

In order to initialize DEXOM, run the following command in the root folder of the project from Matlab:

```
>> dexomInit;
DEXOM: Diversity-based enumeration of optimal context-specific metabolic networks (v0.1.0)
This library was tested using Matlab 2015b, CPLEX v12.8 and COBRA Toolbox v3.0.6.

> Initializing DEXOM library for Matlab
 + Adding external dependencies...
	> Trying to initialize COBRA Toolbox ... Done (already loaded)
> IBM CPLEX selected as the default solver (v128)
> Testing DEXOM ... Done.
> DEXOM is ready to use.
```

During the initialization, DEXOM launches a few quick tests for network reconstruction and enumeration. After getting the previous message, the library is ready to use.

