import os
import argparse
import subprocess
from Bio import SeqIO
from shutil import copyfile
from itertools import repeat
from multiprocessing import Pool, freeze_support


#Run Spades on a directory
def RunSpadesDirectory(inputDir, ouputDir):
    R1List = []
    R2List = []
    outFileList = []
    SpadesFileList = []
    SpadesOutList = []
    BandageOutList = []
    for subdir, dirs, files in os.walk(inputDir):
        R1 = ""
        R2 = ""
        outputFilePath = ""
        SpadesFilePath = ""
        SpadesOutDir = ""
        for file in files:
            if file.endswith(r1_end):
                R1 = os.path.join(subdir, file)
                R2 = os.path.join(subdir, file[:-len(r1_end)]+r2_end)
                R1List.append(R1)
                R2List.append(R2)
                sampleStr = file.replace(r1_end, "")
                outputFilePath = os.path.join(ouputDir, sampleStr)
                outFileList.append(outputFilePath)
                
                SpadesOutDir = os.path.join(ouputDir, sampleStr, "assemble")
                SpadesFilePath = os.path.join(SpadesOutDir, "scaffolds.fasta")
                
                SpadesOutList.append(SpadesOutDir)
                SpadesFileList.append(SpadesFilePath)
                BandageOutList.append(os.path.join(ouputDir, sampleStr, "preview.png"))
    #make out dir for every run
    os.makedirs(os.path.join(ouputDir, sampleStr), 0o777, True)

    RunSpadesParallel(R1List, R2List, SpadesOutList, jobs, threads)
    RunBandageParallel(outFileList, BandageOutList)
    #RunProkkaParallel(SpadesFilePath, outFileList, SpadesFilePath) #prefix?

#Run on outDir's Spades assemble out put    
def AnnotatePhage(Dir):
    prefixList = []
    contigsFileList = []
    contigsOutDirList = []
    for run in os.listdir(Dir):
        for assemble in os.listdir(os.path.join(Dir, run)):
            contigsOutDir = ""
            contigsOutPath = ""
            SpadesFilePath = os.path.join(Dir, run, assemble, "scaffolds.fasta")
            print(SpadesFilePath)
            if os.path.exists(SpadesFilePath):
                for contigs in SeqIO.parse(SpadesFilePath, "fasta"):
                        if len(contigs) > dlen:
                            contigsOutDir = os.path.join(Dir, run ,contigs.name) #OutDir is the sub dir of run
                            contigsOutPath = os.path.join(contigsOutDir, contigs.name+".fasta")
                            contigsOutDirList.append(contigsOutDir)
                            contigsFileList.append(contigsOutPath)
                            prefixList.append(contigs.name)
                            os.makedirs(contigsOutDir, 0o777, True)
                            SeqIO.write(contigs, contigsOutPath, "fasta")
        RunProkkaParallel(contigsFileList, contigsOutDirList, prefixList)

def RunDiamondDir(Dir):
    filePathList = []
    outFileList = []
    for subdir, dirs, files in os.walk(Dir):
        filePath = ""
        for file in files:
            if file.endswith("faa") and "NODE" in file:
                filePath = os.path.join(subdir, file)
                outFile = os.path.join(subdir, file, ".tsv")
                filePathList.append(filePath)
                outFileList.append(outFile)
    RunDiamondParallel(filePathList, outFileList)



#Run Spades in parallel
def RunSpadesParallel(R1List, R2List, outFileList, jobs, threads):
    #numOfprocess = len(R1List)
    #pool = Pool(processes=numOfprocess)
    pool = Pool(processes=jobs)
    pool.starmap(RunSpades, zip(R1List, R2List, outFileList, repeat(threads)))
    pool.close()
    pool.join()
    pool.terminate()

#SPAdes Assembling
def RunSpades(R1, R2, OutDir, threads):
    os.makedirs(OutDir, 0o777, True)
    #cmd = "spades.py --isolate -1 " + R1 + " -2 " + R2 + " -o " + OutDir
    cmd = "spades.py --meta -1 " + R1 + " -2 " + R2 + " -o " + OutDir + " -t " + str(threads)
    subprocess.call(cmd, shell=True)


#Run Bandage in parallel
def RunBandageParallel(fileList, outFileList):
    #numOfprocess = len(R1List)
    #pool = Pool(processes=numOfprocess)
    pool = Pool(processes=4)
    pool.starmap(RunBandage, zip(fileList, outFileList))
    pool.close()
    pool.join()
    pool.terminate()
#Bandage image CD1382_FDSW202399938-1r/assembly_graph.fastg CD1382.jpg
#Bandage Preview
def RunBandage(InFile, OutFile):
    #cmd = "spades.py --isolate -1 " + R1 + " -2 " + R2 + " -o " + OutDir
    cmd = "Bandage image " + InFile + "/assembly_graph.fastg "+ OutFile
    subprocess.call(cmd, shell=True)

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
    cmd = "diamond blastp --db /home/malab/databases_of_malab/nr/nr --query " + fasta + " --out " + outFile + " --evalue 1e-05 --outfmt 6 --max-target-seqs 1 --threads 10"
    subprocess.call(cmd, shell=True)


parser = argparse.ArgumentParser(description='Assembly reads')
parser.add_argument('-i', '--input', dest='fileDir', type=str, required=True,
                    help="the path of the reads")
parser.add_argument('-o', '--output', dest='OpDir', type=str, required=True,
                    help="the output path of reads")
parser.add_argument('-j', '--jobs', dest='jobs', type=str,  required=False, default='4',
                    help="the number of jobs run in parallel")
parser.add_argument('-t', '--threads', dest='threads', type=str, required=False, default='6',
                    help="the number of threads run for a job")
parser.add_argument('-l', '--length', dest='length', type=str, required=False, default='10000',
                    help="the length to filter contigs")
parser.add_argument('-F', '--sepF', dest='sp1', type=str, required=False, default='_1.clean.fq.gz',
                    help="It is the surfix to recognize the forward info, default='_1.clean.fq.gz'.")
parser.add_argument('-R', '--sepR', dest='sp2', type=str, required=False, default='_2.clean.fq.gz',
                    help="It is the surfix to recognize the reverse info, default='_2.clean.fq.gz'.")
args = parser.parse_args()

inputDir = str(args.fileDir)
ouputDir = os.path.abspath(args.OpDir)
jobs = int(args.jobs)
threads = int(args.threads)
#definate length of a phage genome
dlen = int(args.length)
r1_end = args.sp1
r2_end = args.sp2

RunSpadesDirectory(inputDir, ouputDir)
AnnotatePhage(ouputDir)
RunDiamondDir(ouputDir)