# popshift-ms-data
A repository containing data for the initial PopShift manuscript with T4 Lysozyme L99A.

More specifically, it contains the intermediate files and plotting scripts necessary to reproduce the figures from [PopShift: A Thermodynamically Sound Approach to Estimate Binding Free Energies by Accounting for Ligand-Induced Population Shifts from a Ligand-Free Markov State Model](10.1021/acs.jctc.3c00870).

## What does that manuscript describe?

The manuscript describes an alternative approach to assessing binding affinities of a small molecule to an ensemble. It's an approach that allows one to use post-hoc affinity estimations to the states of an MSM in lieu of running independent free energy or docking calculations to complete simulations.

## Where can I get help with this method?

As with all free and open source academic software, I make no guarantee of support for these scripts, or aid for your project. If you are having issues applying PopShift to your problem, or especially if you have found a bug in PopShift, please raise an [issue](https://github.com/bowman-lab/PopShift/issues) in that repository. If you are looking for help understanding what's going on here, please instead raise an [issue](https://github.com/bowman-lab/popshift-ms-data/issues) in this repository.

While it's always hard to reconstruct what one has done over the course of more than a year, I believe that what has been transferred here has not been broken in the process of transfer. If you do try to use some of what I've done here, either for validation or as a suggestion for your own work, and you find something doesn't work it will likely be for two reasons:
1. Path issues: for example, the shell script you're trying has a path to a python script that is different to the one I used on my cluster, and so it needs to be modified to match yours.
2. A file is missing: for example, my script reads a text file that didn't get copied into the repository.

I've made an effort to guard against the second type of issue, but because this isn't software per se and I don't know how things will be set up on your system the first issue is much harder for me to address. What's more, if one knows that issues of the first kind may emerge it's easier to troubleshoot them. Either type of problem would be interesting to hear about if you feel that you are at an impasse--there may be some way for me to help.

## What is the basic workflow these scripts perform?

You should check out the github page for [PopShift](https://github.com/bowman-lab/popshift) to see a general description of the workflow. In order to run some of the scripts in this repository you will need to clone that repository and install its dependencies (easy to do with conda!) in any case. However, I will also mention a few things that are probably at best described as expansions of the workflow discussed there, since they were needed for the paper. If you are here primarily to see whether popshift can be helpful to you, I recommend you not use these alternative procedures since many are excessive, computationally expensive, and unneeded. If you are trying to validate what I did, or are trying to validate some new variant of the method, they might be helpful places to start. The following subsections describe each task the scripts perform, in the order they need to be performed.

### Make Markov state models from simulation datasets

- Run simulations of lysozyme, and prepare three sets of these for msm construction.[^1]
- Three sets of simulation data were generated completely independently, so results associated with them are saved in directories labeled `t4l-X`, where `X` is the identity of the dataset.
- Using `dihedral_tica_featurize.py` and `featurize.sbatch`, obtain dihedral angle featurizations of each trajectory (note in my case these are made with pocket dihedrals only).
- Perform tica and reduce dimensions with `dihedral_tica_reduce.py` and `tica_reduce.sbatch`
- use 'VAMP-2' scan to determine the lighest number of clusters that isn't overfitted with `parallel_vamp_scan_train_test_split.py`, and example usage in 'vamp.sbatch'.
- Plot the results with `plot-vamps-replicas-pocket.py`
- Pick a number of clusters based on a divergence between train and test data in these plots--in our case 75 was chosen.
- perform final clustering with full dataset with that number of clusters, and `dihedral_tica_clustering.py` and `clustering.sbatch`
- Fit a model using implied timescales as a guide with `build-msms.py`, saving models to `t4l-X/models`.

### Extract frames from a simulation dataset with accompanying MSM

- Use `pick_align_frames.py` to rifle through each simulation dataset and extract frames for docking.
- Each dataset was truncated by multiple resect factors (see Fig 3 from paper)--0.2, 0.4, 0.6, 0.8, and 1.0. Each truncation was done by taking the first y frames of the trajectory, where y is trajectory length times the resect factor. Frames were extracted independently at random.
- each extraction was organized by state index from clustering the aforementioned trajectories.
- Multiple extrations were performed on some resects to test the impact of extracting different numbers of frames per state.
- These choices are recorded in the file paths for the extracted frames. For example `t4l-3/binding/resect-1.0/tica-500-msm-2000-k-75-nframes-20/receptor/0/1-332501.pdb` is the path to a frame extracted:
  - from replica 3
  - from resect-1.0, a clustering made with the full-length dataset
  - from a clustering made on a tica decomposition (retaining 90% of kinetic variance) with a 500 frame, 5 ns lagtime, a 2000 frame or 20 ns lagtime msm, with 75 k-means states.
  - from frames categorized as belonging to state `0`.
  - trajectory number 1, frame 332501.

### Prepare these files for docking

- script uses gnu [`parallel`](https://www.gnu.org/software/parallel/) to orchestrate calling `prepare_receptor.py` from the autodock vina conda package on all receptor files.

### Dock to these frames using [smina](https://github.com/mwojcikowski/smina)

- The orchestration of this multitude of docking jobs is done by `dock-t4l-resects.sbatch`
- docked poses can be found in a directory structure that matches that of the cognate conformation in `receptor`.
  - The cognate pose to the receptor conformation discussed [above](#extract-frames-from-a-simulation-dataset-with-accompanying-msm) would be: `t4l-3/binding/resect-1.0/tica-500-msm-2000-k-75-nframes-20/12xsmina/0/1-332501.pdb`
  - Note that `12xsmina` is a nickname given to the docking run used in the paper because it was performed with a 12 angstrom box centered at 0,0,0 with smina.
  - Because the alignment process pick-

### Process the docking scores

- Extract scores and perform popshift calculation with `extract-and-calx-resects.sh`
- This script calls `extract-scores.py` to collate the scores from each state as an `enspara` ragged array.
- This ragged array is then analyzed using `popshift.py` to obtain extracted score jsons. 
- The paths of each of these follow the pattern described above for the receptor structures.

### Create plots for figures

[^1]: Note that the trajectories are not here because they'd be too large, but I'm trying to find a way to host them.
