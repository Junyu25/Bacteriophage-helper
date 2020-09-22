import os
import argparse
import subprocess
from Bio import SeqIO
from shutil import copyfile
from itertools import repeat
from multiprocessing import Pool, freeze_support

#Run on Dir    
def RunProkkaDir(Dir):
    prefixList = []
    fileList = []
    outFileList = []
    for run in os.listdir(Dir):
        if os.path.exists(os.path.join(Dir, run)):
            prefixList.append(run.replace(".fasta", ""))
            fileList.append(os.path.join(Dir, run))
            outFileList.append(os.path.join(ouputDir, run.replace(".fasta", "")))
    RunProkkaParallel(fileList, outFileList, prefixList)

def RunDiamondDir(Dir):
    filePathList = []
    outFileList = []
    for subdir, dirs, files in os.walk(Dir):
        filePath = ""
        for file in files:
            if file.endswith("faa") and "NODE" in file:
                filePath = os.path.join(subdir, file)
                outFile = os.path.join(subdir, file + ".tsv")
                filePathList.append(filePath)
                outFileList.append(outFile)
    RunDiamondParallel(filePathList, outFileList)



#Run Prokka in parallel
def RunProkkaParallel(fileList, outFileList, prefixList):
    #numOfprocess = len(R1List)
    #pool = Pool(processes=numOfprocess)
    pool = Pool(processes=4)
    pool.starmap(RunProkka, zip(fileList, outFileList, prefixList))
    pool.close()
    pool.join()
    pool.terminate()

def RunProkka(fasta, outDir, prefix):
    cmd = "prokka --kingdom Viruses --genus viral --hmms /home/junyuchen/Lab/Phage-SOP/Database/VOGDB/VOGDB_m.hmm --addgenes --prefix " + prefix + " --outdir " + outDir + " --force " + fasta
    print(cmd)
    subprocess.call(cmd, shell=True)

#Run Diamond in parallel
def RunDiamondParallel(fileList, outFileList):
    #numOfprocess = len(R1List)
    #pool = Pool(processes=numOfprocess)
    pool = Pool(processes=2)
    pool.starmap(RunDiamond, zip(fileList, outFileList))
    pool.close()
    pool.join()
    pool.terminate()

def RunDiamond(fasta, outFile):
    cmd = "diamond blastp --db " + nr + " --query " + fasta + " --out " + outFile + " --evalue 1e-05 --outfmt 6 --max-target-seqs 1 --threads 10"
    subprocess.call(cmd, shell=True)


parser = argparse.ArgumentParser(description='Phage Assembly & Annotation')
parser.add_argument('-i', '--input', dest='fileDir', type=str, required=True,
                    help="the path of the reads")
parser.add_argument('-o', '--output', dest='OpDir', type=str, required=True,
                    help="the output path of reads")
parser.add_argument('-d', '--database', dest='Database', type=str, required=False, default='/home/malab/databases_of_malab/nr/nr',
                    help="the nr database path")
parser.add_argument('-j', '--jobs', dest='jobs', type=str,  required=False, default='4',
                    help="the number of jobs run in parallel")
parser.add_argument('-t', '--threads', dest='threads', type=str, required=False, default='6',
                    help="the number of threads run for a job")
args = parser.parse_args()

inputDir = os.path.abspath(args.fileDir)
ouputDir = os.path.abspath(args.OpDir)
nr = os.path.abspath(args.Database)
jobs = int(args.jobs)
threads = int(args.threads)

RunProkkaDir(inputDir)