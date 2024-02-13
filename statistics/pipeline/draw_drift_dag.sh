export PYTHONPATH=~/.venv/lib64/python3.7/site-packages

## activate virtual environment
~/.venv/bin/python3 ~/.venv/bin/activate_this.py 

~/.venv/bin/snakemake --snakefile /home/byrne/halo/dev/hodgkins_dev/pipeline/Drift_Snakefile -n --forceall --rulegraph run_drift | dot -Tsvg > hodgkins_drift_rules.svg
~/.venv/bin/snakemake --snakefile /home/byrne/halo/dev/hodgkins_dev/pipeline/Drift_Snakefile -n --forceall --dag run_drift | dot -Tsvg > hodgkins_drift_dag.svg

