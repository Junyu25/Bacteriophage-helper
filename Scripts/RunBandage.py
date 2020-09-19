import os
import sys
import argparse
import subprocess
#from Bio import SeqIO
from shutil import copyfile
from itertools import repeat
from multiprocessing import Pool, freeze_support


#Run Spades on a directory
def RunBandageDir(inputDir):
    SpadesFileList = []
    BandageOutList = []
    for Dir in os.listdir(inputDir):
        SpadesFilePath = ""
        if os.path.exists(os.path.join(inputDir, Dir, "assemble", "assembly_graph.fastg")):
            SpadesFilePath = os.path.join(inputDir, Dir, "assemble", "assembly_graph.fastg")
            SpadesFileList.append(SpadesFilePath)
            #Bandage
            BandageOutList.append(os.path.join(inputDir, Dir, Dir+"_preview.png"))

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
    cmd = "Bandage image " + InFile + " "+ OutFile
    subprocess.call(cmd, shell=True)



inputDir = sys.argv[1]
RunBandageDir(inputDir)