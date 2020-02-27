1 #!/bin/bash
      2 #SBATCH -p short
      3 #SBATCH --open-mode=truncate
      4 #SBATCH --mail-type=FAIL
      5 #SBATCH --mail-user=alanlegoallec@g.harvard.edu
      6
      7 module load gcc/6.2.0
      8 module load python/3.6.