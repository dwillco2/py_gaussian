Various scripts for the automation of quantum chemical jobs using gaussian, orca, and CREST.

Each job type has its own class, e.g. the GenGaussian class which contains various utilities for automated job submission and inherit from a GenJob class.

Currently set up for jobs on the Heriot-Watt DMOG cluster using slurm. Can be easily modified to suit your HPC environment.
