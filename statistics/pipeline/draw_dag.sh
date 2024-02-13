#### MAKE SURE SNAKEFILE DOES NOT PRINT ANY MESSAGES, AS THEY WILL CAUSE ERRORS 
#### FOR DOT COMMAND

export PYTHONPATH=~/.venv/lib64/python3.7/site-packages

## activate virtual environment
~/.venv/bin/python3 ~/.venv/bin/activate_this.py 

~/.venv/bin/snakemake --snakefile /home/byrne/halo/dev/hodgkins_dev/pipeline/Snakefile -n --forceall --rulegraph all | dot -Tsvg > hodgkins_rules.svg
#~/.venv/bin/snakemake --snakefile /home/byrne/halo/dev/hodgkins_dev/pipeline/Snakefile -n --forceall --dag all | dot -Tsvg > hodgkins_dag.svg

