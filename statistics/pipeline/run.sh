export PYTHONPATH=/home/byrne/.venv/lib64/python3.7/site-packages
export LSB_JOB_REPORT_MAIL=Y

## activate virtual environment
/home/byrne/.venv/bin/python3 /home/byrne/.venv/bin/activate_this.py 

bsub -J 'hodgkins_pipeline' -e 'hodgkins_pipeline.log' -o 'hodgkins_pipeline.log' \
  -W 3000 \
  ~/.venv/bin/snakemake --jobs 499 --cluster '/home/byrne/halo/dev/hodgkins_dev/pipeline/bsub.py' \
    --configfile input/config/study_config.yaml \
    --cluster-config /home/byrne/halo/dev/hodgkins_dev/pipeline/lsf.yaml \
    --snakefile /home/byrne/halo/dev/hodgkins_dev/pipeline/Snakefile \
    -p -r "$@" 

