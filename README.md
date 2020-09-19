# Bacteriophage-helper
Bacteriophage standard operating procedure

## Usage

First Activate the Environment:
`source /home/junyuchen/Biosoft/anaconda3/bin/activate /home/junyuchen/Biosoft/anaconda3/envs/phage`

Second Use the pipeline:
```
usage: phage-sop.py [-h] -i FILEDIR -o OPDIR [-j JOBS] [-t THREADS] [-l LENGTH] [-F SP1] [-R SP2]

Phage Assembly & Annotation

optional arguments:
  -h, --help            show this help message and exit
  -i FILEDIR, --input FILEDIR
                        the path of the reads
  -o OPDIR, --output OPDIR
                        the output path of reads
  -j JOBS, --jobs JOBS  the number of jobs run in parallel
  -t THREADS, --threads THREADS
                        the number of threads run for a job
  -l LENGTH, --length LENGTH
                        the length to filter contigs
  -F SP1, --sepF SP1    It is the surfix to recognize the forward info, default='_1.clean.fq.gz'.
  -R SP2, --sepR SP2    It is the surfix to recognize the reverse info, default='_2.clean.fq.gz'.
```


eg.
```shell
python /home/junyuchen/Lab/Phage-SOP/Bacteriophage-helper/Scripts/phage-sop.py -i /home/junyuchen/Lab/Phage-SOP/rawdata/jtshen-2020-06-13 -o /home/junyuchen/Lab/Phage-SOP/Out-test
```
