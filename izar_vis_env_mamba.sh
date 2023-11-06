# script to install environment
mamba create --name izar_vis python=3.9
mamba activate izar_vis
mamba install biopython matplotlib numpy pandas seaborn tqdm -c conda-forge -c anaconda
pip install fair-esm