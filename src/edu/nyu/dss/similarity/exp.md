1. copy datasets to scratch

scp -r /Users/sw160/Desktop/spadas-dataset sw160@greene.hpc.nyu.edu:/scratch/sw160/spadas/dataset

scp -r /Users/sw160/Desktop/argoverse-api sw160@greene.hpc.nyu.edu:/scratch/sw160/spadas/dataset

2. update jar and slurm

mvn clean package

./update_greene.sh

3. run with slurm
sbatch runSpadas.s

4. collect logs

scp -r sw160@greene.hpc.nyu.edu:/scratch/sw160/spadas/logs /Users/sw160/Desktop/