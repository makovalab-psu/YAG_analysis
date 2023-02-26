#!/bin/bash
#SBATCH --job-name=marta
#SBATCH --output=marta-%j.out
#SBATCH --error=marta-%j.err
#SBATCH -C new
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=500G


/nfs/brubeck.bx.psu.edu/scratch5/marta/paml4.8/bin/codeml *.ctl -t 64





