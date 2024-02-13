export PYTHONPATH=~/.venv/lib64/python3.7/site-packages

## after reinstalling ImageMagick into my own directory, these were the only two things that needed to be set in the environment.
## did not work with Zheng's version because his was configured and installed in /ifs originally
export PATH=/home/byrne/halo/dev/hodgkins_dev/bin/ImageMagick-7.0.11-7/bin:/home/byrne/halo/dev/hodgkins_dev/drift_detection/opencv-lib-3.4.7-juno/bin:$PATH

export LD_LIBRARY_PATH=/home/byrne/halo/dev/hodgkins_dev/bin/ImageMagick-7.0.11-7/lib:/home/byrne/halo/dev/hodgkins_dev/drift_detection/opencv-lib-3.4.7-juno/lib64:$LD_LIBRARY_PATH

## activate virtual environment
~/.venv/bin/python3 ~/.venv/bin/activate_this.py 

bsub -J 'hodgkins_drift' -e 'hodgkins_drift.log' -o 'hodgkins_drift.log' -W 36000 \
   ~/.venv/bin/snakemake --jobs 299 --cluster '/home/byrne/halo/dev/hodgkins_dev/pipeline/bsub.py' \
    --configfile input/config/study_config.yaml \
    --cluster-config /home/byrne/halo/dev/hodgkins_dev/pipeline/lsf.yaml \
    --snakefile /home/byrne/halo/dev/hodgkins_dev/pipeline/Drift_Snakefile \
    -p -r "$@" > hodgkins_drift.log 2>&1 & #\
#    run_drift
