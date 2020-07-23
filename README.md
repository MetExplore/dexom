# DEXOM 
[![BioRxiv](https://img.shields.io/badge/BioRxiv-2020.07.17.208918-brightgreen)](https://www.biorxiv.org/content/10.1101/2020.07.17.208918v1)
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

## Quick start

The usage of DEXOM is similar to any of the context-specific network reconstruction & data integration methods included in COBRA Toolbox. If you are not familiar with these methods, the tutorial "[Extraction of context-specific models](https://opencobra.github.io/cobratoolbox/stable/tutorials/tutorialExtractionTranscriptomic.html)" is a good starting point.

The method requires at least 3 things to work: a Genome-Scale Metabolic Network (GSMN), a list of reactions associated with highly expressed enzymes, and a list of reactions associated with lowly expressed enzymes. By default, DEXOM uses the same objective function as the one used by the iMAT algorithm, i.e., it tries to extract sub-networks from the provided GSMN maximizing the selection of reactions associated with highly expressed enzymes, and minimizing the inclusion of reactions associated with lowly expressed enzymes.

The library includes some toy models for testing the algorithm. Here is an example to enumerate 5 optimal metabolic sub-networks in the DAG model introduced in the research paper:

```matlab
% Make sure dexom is loaded with the embedded COBRA Toolbox v3.0.6
dexomInit;
% Create a DAG metabolic network with 5 layers and 4 metabolites per layer
model = dagNet(5,4);
enumOptions.maxUniqueSolutions = 5;
% Indexes of the reactions in the model that are associated
% with highly expressed enzymes
methodOptions.RHindex = [];
% Indexes of the reactions in the model that are associated
% with lowly expressed enzymes. For this toy example, all rxn
% are assumed to be associated with lowly abundant enzymes
methodOptions.RLindex = 1:length(model.rxns);
results = dexom(model, methodOptions, enumOptions);
```

The method returns a structure with the results and many other useful output values for posterior analysis of the results. In order to extract the unique solutions as row binary vectors (indicating with 1s the selected reactions from the model), you can use the method `getUniqueAcceptedSolutions`:

```matlab
solutions = getUniqueAcceptedSolutions(results);
% Extract the first optimal context-specific model
selectedRxn = solutions(1,:);
cModel1 = removeRxns(model, model.rxns(selectedRxn == 0));
% Test that the model is flux consistent
cOpts.epsilon=1e-6;
cOpts.modeFlag=0;
cOpts.method='fastcc';
[~,fluxConsistentRxnBool] = findFluxConsistentSubset(cModel1, cOpts);
assert (sum(fluxConsistentRxnBool == 0) == 0)
```

The algorithm is highly customizable. Options are divided in two different structures, one containing the options for the context-specific method, and another with the configuration for the enumeration of alternative optimal solutions. For more information about the possible parameters of the method, see [dexomDefaultOptions.m](https://github.com/MetExplore/dexom/blob/master/src/methods/dexom/dexomDefaultOptions.m), and for the possible parameters to adjust the behavior of the enumeration strategy, see [defaultEnumOptions.m](https://github.com/MetExplore/dexom/blob/master/src/methods/defaultEnumOptions.m)

## Reproducibility of the experiments

The library includes a submodule with the data and all the output files (matlab files and exported csv files). Scripts to reproduce all the steps described in the research paper are also available in the [evaluation]() folder. To clone the repository with all the data used in the experiments as well as the files resulting from the reconstruction, use the following git command:

```
git clone --recursive https://github.com/MetExplore/dexom.git
```

The scripts to reproduce the analysis are:

* [scriptEvaluationSamplingDAG.m](https://github.com/MetExplore/dexom/blob/master/src/evaluation/scriptEvaluationSamplingDAG.m). This script contains the code to sample up to 250 unique solutions in the DAG network model.
* [scriptEvaluationSamplingYeast6.m](https://github.com/MetExplore/dexom/blob/master/src/evaluation/scriptEvaluationSamplingYeast6.m). Contains the code to perform reconstruction using the Yeast 6 model with random sets of highly expressed and lowly expressed genes.
* [scriptYeastEvaluation.m](https://github.com/MetExplore/dexom/blob/master/src/evaluation/scriptYeastEvaluation.m). Script to generate and evaluate the ensembles to predict essential genes in yeast using the Yeast 6 model.

Note that there is some amount of randomness involved in the reconstruction process, depending on the solver version/model and configuration. Small variations of the results are expected between simulations.

## How to cite

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



