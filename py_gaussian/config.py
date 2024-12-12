bash_template = r"""#!/bin/bash -l

################# Part-1 Slurm directives ####################
## Environment variables
#SBATCH --export=ALL
## Output and Error Files
#SBATCH -o job-%j.output
#SBATCH -e job-%j.error
## Job name
#SBATCH -J _JOBNAME_
## Run time: "hours:minutes:seconds", "days-hours"
#SBATCH --time=_RUNTIME_:00:00
## Memory limit (in megabytes). Total --mem or amount per cpu --mem-per-cpu
#SBATCH --mem-per-cpu=_MEMORY_
## Processing slots
#SBATCH --nodes=1
#SBATCH --ntasks=_NPROCS_
## Specify partition
#SBATCH -p nodes

################# Part-2 Run Job ####################
"""