# Generating SS Ricci


1. Copy `create_scripts`, `d2R.wl` and `d2R-SS-r0.wl` to the machine where
   the calculation is to be run.
2. Edit the "Parameters" section of `create_scripts` and run. This will generate
   a set of directories, one for each radius.
3. Submit the jobs to the cluster using the generated `sonic.sh` SLURM scripts.
   When complete, there will be a `d2R-data-r0_XXX.m` file.
4. Run `Convert-d2RSS.nb` to generate HDF5 files.
