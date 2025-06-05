# PhyllosphereDynamics
 Models and scripts for simulating population dynamics from metabolic interactions between leaf-associated bacteria

========================

This repository contains the genome-scale metabolic models reported in Pacheco, Ugolini *et al.*, as well as the scripts required to simulate dynamics stemming from metabolic exchange using [COMETS v2.11.3](https://www.runcomets.org/home).

## Models

Located in the 'Models' directory, contains genome-scale models for *Sphingomonas* sp. Leaf257 and *Rhizobium* sp. Leaf68, two bacterial members of the *Arabidopsis thaliana* phyllosphere microbiome (*At*-LSPHERE, Bai *et al.*, 2015) originally constructed in Schäfer, Pacheco, *et al.*, 2023. The models are provided in .sbml format for use with [COBRApy](https://opencobra.github.io/cobrapy/).

## Simulation scripts and output

Located in the 'Simulation' directory, contains files for simulating dynamics in COMETS between the two models in a minimal medium containing xylan. Time-resolved biomass abundances and fluxes for each model are contained in the 'biomassLog,' and 'fluxLog.m' files, respectively.

The script 'readCometsLogs' can be used to plot the relevant biomass, media, and flux logs. 

## Required dependencies
  * [MATLAB](https://www.mathworks.com/products/matlab.html) R2021a or higher
  * [COBRA Toolbox](https://opencobra.github.io/cobratoolbox/stable/) v2.24.3 or higher
  * [COMETS MATLAB toolbox](https://github.com/segrelab/comets-toolbox)