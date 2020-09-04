import os
import argparse
import subprocess
from Bio import SeqIO
from shutil import copyfile
from itertools import repeat
from multiprocessing import Pool, freeze_support

#definate length of a phage genome
dlen = 10000

#Run Spades on a directory
def RunSpadesDirectory(inputDir, ouputDir):
    for subdir, dirs, files in os.walk(inputDir):
        R1 = ""
        R2 = ""
        outputFilePath = ""
        SpadesFilePath = ""
        SpadesOutDir = ""
        R1List = []
        R2List = []
        outFileList = []
        SpadesFileList = []
        SpadesOutList = []
        BandageOutList = []
        for file in files:
            if file.endswith("_1.clean.fq.gz"):
                R1 = os.path.join(subdir, file)
                R2 = os.path.join(subdir, file[:-13]+"2.clean.fq.gz")
                R1List.append(R1)
                #print(R1)
                #print(R2)
                R2List.append(R2)
            sampleStr = os.path.splitext(file)[0]
            outputFilePath = os.path.join(ouputDir, sampleStr)
            outFileList.append(outputFilePath)
            
            SpadesOutDir = os.path.join(ouputDir, sampleStr, "assemble")
            SpadesFilePath = os.path.join(SpadesOutDir, "scaffolds.fasta")
            
            SpadesOutList.append(SpadesOutDir)
            SpadesFileList.append(SpadesFilePath)
            BandageOutList.append(os.path.join(ouputDir, sampleStr, "preview.png"))
    #make out dir for every run
    os.makedirs(os.path.join(ouputDir, sampleStr), 0o777, True)

    RunSpadesParallel(R1List, R2List, SpadesOutList)

    RunBandageParallel(outFileList, BandageOutList)
    
    RunProkkaParallel(SpadesFilePath, outFileList, SpadesFilePath) #prefix?
#Run on outDir's Spades assemble out put    
def AnnotatePhage(Dir):
    for run in os.listdir(Dir):
        for assemble in os.listdir(os.path.join(Dir, run)):
            contigsOutDir = ""
            contigsOutPath = ""
            #SpadesFilePath = ""
            prefixList = []
            contigsFileList = []
            contigsOutDirList = []
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
                            #print(contigsOutDir)
                            #print(contigsOutPath)
                            os.makedirs(contigsOutDir, 0o777, True)
                            SeqIO.write(contigs, contigsOutPath, "fasta")
        RunProkkaParallel(contigsFileList, contigsOutDirList, prefixList)
        

'''      
def PreContigs(contigs, OutDir):
    for contigs in SeqIO.parse(contigs, "fasta"):
    if len(contigs) > dlen:
        print(contigs.id)
        print(contigs)
        contigsOutDir = os.path.join(OutDir, contigs.name)
        contigsOutPath = os.path.join(contigsOutDir, contigs.name+".fasta")
        os.makedirs(contigsOutDir, 0o777, True)
        SeqIO.write(contigs, contigsOutPath, "fasta")
'''        
'''        
def CopyResult(InList, OutList):
    #copy assemble result to parent dir
    copyfile(os.path.join(SpadesOutDir, "scaffolds.fasta"), outputFilePath)
'''

#Run Spades in parallel
def RunSpadesParallel(R1List, R2List, outFileList):
    #numOfprocess = len(R1List)
    #pool = Pool(processes=numOfprocess)
    pool = Pool(processes=20)
    pool.starmap(RunSpades, zip(R1List, R2List, outFileList))
    pool.close()
    pool.join()
    pool.terminate()

#SPAdes Assembling
def RunSpades(R1, R2, OutDir):
    os.makedirs(OutDir, 0o777, True)
    #cmd = "spades.py --isolate -1 " + R1 + " -2 " + R2 + " -o " + OutDir
    cmd = "spades.py --meta -1 " + R1 + " -2 " + R2 + " -o " + OutDir
   
    subprocess.call(cmd, shell=True)


#Run Bandage in parallel
def RunBandageParallel(fileList, outFileList):
    #numOfprocess = len(R1List)
    #pool = Pool(processes=numOfprocess)
    pool = Pool(processes=20)
    pool.starmap(RunBandage, zip(fileList, outFileList))
    pool.close()
    pool.join()
    pool.terminate()
#Bandage image CD1382_FDSW202399938-1r/assembly_graph.fastg CD1382.jpg
#Bandage Assembling
def RunBandage(InFile, OutFile):
    #cmd = "spades.py --isolate -1 " + R1 + " -2 " + R2 + " -o " + OutDir
    cmd = "Bandage image " + InFile + "/assembly_graph.fastg "+ OutFile
    subprocess.call(cmd, shell=True)

def RunProkkaParallel(fileList, outFileList, prefixList):
    #numOfprocess = len(R1List)
    #pool = Pool(processes=numOfprocess)
    pool = Pool(processes=20)
    pool.starmap(RunProkka, zip(fileList, outFileList, prefixList))
    pool.close()
    pool.join()
    pool.terminate()

def RunProkka(fasta, outDir, prefix):
    cmd = "prokka --kingdom Viruses --genus viral --hmms /home/junyuchen/Lab/Phage-SOP/Database/VOGDB/VOGDB_m.hmm --addgenes --prefix " + prefix + " --outdir " + outDir + " --force " + fasta
    print(cmd)
    subprocess.call(cmd, shell=True)


parser = argparse.ArgumentParser(description='Assembly reads')
parser.add_argument('-i', '--input', dest='fileDir', type=str, required=True,
                    help="the path of the reads")
parser.add_argument('-o', '--output', dest='OpDir', type=str, required=True,
                    help="the output path of reads")
args = parser.parse_args()

inputDir = str(args.fileDir)
ouputDir = os.path.abspath(args.OpDir)

RunSpadesDirectory(inputDir, ouputDir)
AnnotatePhage(ouputDir)