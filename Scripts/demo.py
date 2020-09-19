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
            # print(R1)
            # print(R2)
            R2List.append(R2)
            sampleStr = os.path.splitext(
                file)[0].replace("_1.clean.fq", "")
            outputFilePath = os.path.join(ouputDir, sampleStr)
            outFileList.append(outputFilePath)

            SpadesOutDir = os.path.join(ouputDir, sampleStr, "assemble")
            SpadesFilePath = os.path.join(SpadesOutDir, "scaffolds.fasta")

            SpadesOutList.append(SpadesOutDir)
            SpadesFileList.append(SpadesFilePath)
            BandageOutList.append(os.path.join(
                ouputDir, sampleStr, "preview.png"))
