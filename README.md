# popshift-ms-data
A repository containing data for the initial PopShift manuscript with T4 Lysozyme L99A.

More specifically, it contains the intermediate files and plotting scripts necessary to reproduce the figures from [PopShift: A Thermodynamically Sound Approach to Estimate Binding Free Energies by Accounting for Ligand-Induced Population Shifts from a Ligand-Free Markov State Model](10.1021/acs.jctc.3c00870).

## What does that manuscript describe?

The manuscript describes an alternative approach to assessing binding affinities of a small molecule to an ensemble. It's an approach that allows one to use post-hoc affinity estimations to the states of an MSM in lieu of running independent free energy or docking calculations to complete simulations.

## What is the basic workflow these scripts perform?

You should check out the github page for [PopShift](https://github.com/bowman-lab/popshift) to see a general description of the workflow. In order to run some of the scripts in this repository you will need to clone that repository and install its dependencies (easy to do with conda!) in any case. However, I will also mention a few things that are probably at best described as expansions of the workflow discussed there, since they were needed for the paper. If you are here primarily to see whether popshift can be helpful to you, I recommend you not use these alternative procedures since many are excessive, computationally expensive, and uneeded.

These scripts do the following things:
 - Make MSMs from a simulation dataset
   - Using `dihedral_tica_featurize.py` and `featurize.sbatch`, obtain dihedral angle featurizations of each trajectory.
   - Perform tica and reduce dimensions with `dihedral_tica_reduce.py` and `tica_reduce.sbatch`
   - use 'VAMP-2' scan to determine the lighest number of clusters that isn't overfitted.
 - Extract frames from a simulation dataset that was too large to upload to this repository but which I hope to host elsewhere.
   - That dataset was three sets of simulations performed independently as replicas. Data associated with each replica will be in `t4l-X/` where x is the replica number.
   - Each dataset was truncated by multiple resect factors--0.2, 0.4, 0.6, 0.8, and 1.0. Each truncation was done by taking the first y frames of the trajectory, where y is trajectory length times the resect factor. Frames were extracted independently at random.
   - each extraction was organized by state index from clustering the aforementioned trajectories.
   - Multiple extrations were performed on some resects to test the impact of extracting different numbers of frames per state.
   - These choices are recorded in the file paths for the extracted frames. For example `t4l-3/binding/resect-1.0/tica-500-msm-2000-k-75-nframes-20/receptor/0/1-332501.pdb` is the path to a frame extracted:
     - from replica 3
     - from resect-1.0, a clustering made with the full-length dataset
     - from a clustering made on a tica decomposition (retaining 90% of kinetic variance) with a 500 frame, 5 ns lagtime, a 2000 frame or 20 ns lagtime msm, with 75 k-means states.
     - from frames categorized as belonging to state `0`.
     - trajectory number 1, frame 332501.
- Prepare these files for docking using prep-receptor
- Dock to these frames using [smina](https://github.com/mwojcikowski/smina). The orchestration of this multitude of docking jobs is done by `dock-t4l-resects.sbatch`
- Process the docking scores