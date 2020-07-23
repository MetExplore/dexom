# DEXOM
> Diversity-based enumeration of optimal context-specific metabolic networks

DEXOM is a Matlab library for the reconstruction and enumeration of diverse optimal context-specific metabolic networks. It requires COBRA Toolbox (included as submodule) and a MILP solver (CPLEX, Gurobi).

<p align="center"><img src="https://github.com/MetExplore/dexom/raw/master/assets/overview.png" width="500"></p>

## Installation

In order to use DEXOM, you need a compatible version of Matlab (recommended version >= 2015), CPLEX (recommended version >= 12.8), and COBRA Toolbox (recommended v3.0.6). This repository includes an embedded version of COBRA Toolbox v3.0.6 with minor changes to facilitate the reproducibility of the experiments and avoid automatic updates that can introduce incompatibilities. In order to clone DEXOM with the embedded COBRA v3.0.6, use the following git command:

```
git clone --recurse-submodules="modules/cobratoolbox" https://github.com/MetExplore/dexom.git
```

Once cloned, open the dexom folder in matlab and run the initialization script `dexomInit.m`. By default, it removes from your Matlab's path your current COBRA Toolbox and replaces it with the embedded COBRA Toolbox v3.0.6. This is done in order to prevent incompatibilities with other versions. If you want to try to set up everything using your own COBRA Toolbox installation, run `dexomInit(1)`. You can also automatically replace the embedded version with your previous COBRA Toolbox by running `restoreCobraToolboxPath` after initializing DEXOM having used `dexomInit`.

```
>> dexomInit;

DEXOM: Diversity-based Extraction of Optimal Metabolic-networks (v0.1.0)
This version was tested with Matlab 2015b (CPLEX v12.8), 2018a (CPLEX v12.9), 
2018b (CPLEX v12.8) and COBRA Toolbox v3.0.6 on Windows 10.

> Initializing DEXOM library for Matlab
 + Adding external dependencies...
 + Checking and replacing previous COBRA Toolbox installations... 0 entries removed.
 + Initializing the embedded COBRA Toolbox (use dexomInit(0,1) to show log)...
> IBM CPLEX selected as the default solver (v128)
> Testing DEXOM (solver ibm_cplex) ... Done.
> DEXOM is ready to use.
```

During the initialization, DEXOM launches a few quick tests for network reconstruction and enumeration. After getting the previous message, the library is ready to use.

## Citation

If you find this software useful, please consider citing it as:

> Rodriguez-Mier, P., Poupin, N., de Blasio, C., Le Cam, L. & Jourdan, F. DEXOM: Diversity-based enumeration of optimal context-specific metabolic networks. BioRxiv **[Preprint]**. July 17, 2020. Available from: https://doi.org/10.1101/2020.07.17.208918

```
@article {RodriguezMier2020,
	author = {Rodr{\'\i}guez-Mier, Pablo and Poupin, Nathalie and de Blasio, Carlo and Le Cam, Laurent and Jourdan, Fabien},
	title = {DEXOM: Diversity-based enumeration of optimal context-specific metabolic networks},
	year = {2020},
	doi = {10.1101/2020.07.17.208918},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2020/07/17/2020.07.17.208918},
	journal = {bioRxiv}
}
```



