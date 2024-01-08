# What is popshift-ms-data
A repository containing data for the initial PopShift manuscript with T4 Lysozyme L99A.

More specifically, it contains the intermediate files and plotting scripts necessary to reproduce the figures from [PopShift: A Thermodynamically Sound Approach to Estimate Binding Free Energies by Accounting for Ligand-Induced Population Shifts from a Ligand-Free Markov State Model](10.1021/acs.jctc.3c00870).

## What does that manuscript describe?
The manuscript describes an alternative approach to assessing binding affinities of a small molecule to an ensemble. It's an approach that allows one to use post-hoc affinity estimations to the states of an MSM in lieu of running independent free energy or docking calculations to complete simulations.

## What is the basic workflow these scripts perform?
You should check out the github page for [PopShift](https://github.com/bowman-lab/popshift) to see a general description of the workflow. In order to run some of the scripts in this repository you will need to clone that repository and install its dependencies (easy to do with conda!) in any case. However, I will also mention a few things that are probably at best described as expansions of the workflow discussed there, since they were needed for the paper. If you are here primarily to see whether popshift can be helpful to you, I recommend you not use these alternative procedures since many are excessive, computationally expensive, and uneeded.