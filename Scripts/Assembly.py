import os
import argparse
import subprocess
from itertools import repeat
from multiprocessing import Pool, freeze_support

#Run Spades on a directory
def RunSpadesDirectory(inputDir, ouputDir):
    for subdir, dirs, files in os.walk(inputDir):
        R1 = ""
        R2 = ""
        outputFilePath = ""
        R1List = []
        R2List = []
        outFileList = []
        for file in files:
            if file.endswith("_1.clean.fq") and os.path.getsize(file) > 0:
                R1 = os.path.join(subdir, file)
                R1List.append(R1)
            if file.endswith("_2.clean.fq") and os.path.getsize(file) > 0:
                R2 = os.path.join(subdir, file)
                R2List.append(R2)
            sampleStr = os.path.splitext(file)[0]
            outputFilePath = os.path.join(ouputDir, sampleStr)
            outFileList.append(outputFilePath)
    RunSpadesParallel(R1List, R2List, outFileList)

#Run Spades in parallel
def RunSpadesParallel(R1List, R2List, outFileList):
    numOfprocess = len(R1List)
    pool = Pool(processes=numOfprocess)
    #pool = Pool(processes=16)
    pool.starmap(RunSpades, zip(R1List, R2List, outFileList))
    pool.close()
    pool.join()
    pool.terminate()

#SPAdes Assembling
def RunSpades(R1, R2, OutDir):
    os.makedirs(OutDir, 0o777, True)
    cmd = "spades.py -1 " + R1 + " -2 " + R2 + " -o " + OutDir
    subprocess.call(cmd, shell=True)


parser = argparse.ArgumentParser(description='Assembly reads')
parser.add_argument('-i', '--input', dest='fileDir', type=str, required=True,
                    help="the path of the reads")
parser.add_argument('-o', '--output', dest='OpDir', type=str, required=True,
                    help="the output path of reads")
args = parser.parse_args()

inputDir = str(args.file)
ouputDir = os.path.abspath(args.OpDir)

RunSpadesDirectory(inputDir, ouputDir)