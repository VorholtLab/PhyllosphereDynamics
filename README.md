# PhyllosphereDynamics
 Models and scripts for simulating population dynamics from metabolic interactions between leaf-associated bacteria

========================

This repository contains the genome-scale metabolic models reported in Pacheco, Ugolini *et al.*, script for processing simulation logs, and the files for simulating dynamics stemming from metabolic exchange using [COMETS v2.11.3](https://www.runcomets.org/home).

## Models

Located in the 'Models' directory, contains genome-scale models for *Sphingomonas* sp. Leaf257 and *Rhizobium* sp. Leaf68, two bacterial members of the *Arabidopsis thaliana* phyllosphere microbiome (*At*-LSPHERE, Bai *et al.*, 2015) originally constructed in Schäfer, Pacheco, *et al.*, 2023. The models are provided in .sbml format for use with [COBRApy](https://opencobra.github.io/cobrapy/) and in .mat format for use with [the COBRA Toolbox](https://opencobra.github.io/cobratoolbox/stable/index.html). The genome annotations used to generate the corresponding reconstructions in Schäfer, Pacheco, *et al.* are also included in the 'Models/Genomes' directory. 

## Simulation scripts and output

### Input

Located in the 'Simulation' directory, contains files for a representative COMETS simulation between the two models in a minimal medium containing xylan. Following the installation and launching of [COMETS](https://www.runcomets.org/home), this directory is to be set as the location in which to carry out the simulation based on the models and layout files contained within. Please view the [COMETS documentation](https://segrelab.github.io/comets-manual/) for more information on installing and running COMETS.

### Output

The resulting time-resolved biomass abundances, media abundances, and fluxes for each model are contained in the 'biomassLog,' 'mediaLog,' and 'fluxLog.m' files, respectively.

The MATLAB script 'readCometsLogs' is used to plot the relevant biomass, media, and flux logs. 

## System requirements

### Hardware requirements

The plotting script contained within this repository requires only a standard computer with enough RAM to support the in-memory operations. Please view the [COMETS documentation](https://segrelab.github.io/comets-manual/) for more information on hardware requirements for COMETS.

### Required dependencies
  * [MATLAB](https://www.mathworks.com/products/matlab.html) R2021a or higher
  * [COBRA Toolbox](https://opencobra.github.io/cobratoolbox/stable/) v2.24.3 or higher
  * [COMETS MATLAB Toolbox](https://github.com/segrelab/comets-toolbox)
  * [COMETS](https://www.runcomets.org/home) is required to run the simulation within the 'Simulation' directory.

### OS requirements

This script has been tested on: macOS Sequoia (15.5)
The expected runtime of the plotting script is ~2 minutes for a standard desktop computer.