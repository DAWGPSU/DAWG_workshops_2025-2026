Data source: PRJNA880162 ()
Please download samples and place to   `workshop_4/data/PRJNA880162/input_data` folder.

The dataset was processed acording to the instructions, given in workshop_1. Please, follow the instructions of the workshop_1 to download UNITE database and place it to `workshop_4/data/database` directory. To launch the processing pipelines, please use: 
snakemake --snakefile Snakemodule \
          --configfile config_16S.yaml \
          --cluster-config cluster_config.yaml \
          --cluster "sbatch -o {cluster.log} --partition=open --account=open --mem={cluster.mem} --time={cluster.time} --nodes={cluster.nodes} --ntasks={cluster.ntasks}" \
          --use-conda \
          --conda-frontend conda \
          --conda-prefix ../../../conda_envs \
          -p \
          --jobs 60 \
          --cores 10 \
          --latency-wait 60 

snakemake --snakefile Snakemodule \
          --configfile config_ITS.yaml \
          --cluster-config cluster_config.yaml \
          --cluster "sbatch -o {cluster.log} --partition=open --account=open --mem={cluster.mem} --time={cluster.time} --nodes={cluster.nodes} --ntasks={cluster.ntasks}" \
          --use-conda \
          --conda-frontend conda \
          --conda-prefix ../../../conda_envs \
          -p \
          --jobs 60 \
          --cores 10 \
          --latency-wait 60 