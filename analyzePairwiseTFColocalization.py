import numpy as np
import copy
import copy,os,sys,time,datetime,gzip
from tqdm.notebook import tqdm
import multiprocessing

import pandas


def infoLine(message,infoType="info"):
    infoType = infoType.upper()
    if len(infoType) < 5:
        infoType=infoType + " "
    time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    outline = "[" + infoType + " " + str(time) + "] " + message
    print(outline)

    if infoType == "ERROR":
        sys.exit()
#


def processData(infile, tabCountFile, tabFractFile):
    summitDataHash = {}
    factorNameHash = {}
    summitNumByFactor = {}
    cooCountHash = {}
    cooRatioHash = {}
    
    with gzip.open(infile, "rt") as fi:
        for line in fi:
            line = line.rstrip()
            row = line.split("\t")
            
            if len(row) < 8:
                continue
                
            summitName = row[3]
            factorName = row[7]
            
            if summitName not in summitDataHash:
                summitDataHash[summitName] = {}
            summitDataHash[summitName][factorName] = ""
            
            factorNameHash[factorName] = ""
    
    # calculate pairwise cooccurance count
    for summitName in summitDataHash:
        factorNameList = list(summitDataHash[summitName].keys())
        
        for primaryFactor in factorNameList:
            # count factor
            if primaryFactor not in summitNumByFactor:
                summitNumByFactor[primaryFactor] = 1
            else:
                summitNumByFactor[primaryFactor] = summitNumByFactor[primaryFactor] + 1
            
            # count co-occurance
            if primaryFactor not in cooCountHash:
                cooCountHash[primaryFactor] = {}
            
            for affiliatedFactor in factorNameList:
                if affiliatedFactor not in cooCountHash[primaryFactor]:
                    cooCountHash[primaryFactor][affiliatedFactor] = 1
                else:
                    cooCountHash[primaryFactor][affiliatedFactor] = cooCountHash[primaryFactor][affiliatedFactor]  + 1
    # convert count to fraction/ratio
    for primaryFactor in factorNameHash:
        totalCount = 1.0 * summitNumByFactor[primaryFactor]
        cooRatioHash[primaryFactor] = {}
        for affiliatedFactor in cooCountHash[primaryFactor]:
            cooRatioHash[primaryFactor][affiliatedFactor] = int(cooCountHash[primaryFactor][affiliatedFactor] * 10000.0 / totalCount)/10000.0
    
    # output fraction data as table
    factorNameList = list(factorNameHash.keys())
    factorNameList.sort()
    with open(tabFractFile, "wt") as fo:
        fo.write( "primaryFactor\taffiliatedFactor\tpercentage\n" )
        for primaryFactor in factorNameList:
            for affiliatedFactor in factorNameList:
                if primaryFactor == affiliatedFactor:
                    fo.write( primaryFactor + "\t" + affiliatedFactor + "\t" + "1.0" + "\n" )
                else:
                    if affiliatedFactor in cooRatioHash[primaryFactor]:
                        fo.write( primaryFactor + "\t" + affiliatedFactor + "\t" + str(cooRatioHash[primaryFactor][affiliatedFactor]) + "\n" )
                    else:
                        fo.write( primaryFactor + "\t" + affiliatedFactor + "\t" + "0.0" + "\n" )
    # done
    # output count data as table
    factorNameList = list(factorNameHash.keys())
    factorNameList.sort()
    with open(tabCountFile, "wt") as fo:
        fo.write( "primaryFactor\taffiliatedFactor\tcount\n" )
        for primaryFactor in factorNameList:
            for affiliatedFactor in factorNameList:
                if affiliatedFactor in cooCountHash[primaryFactor]:
                    fo.write( primaryFactor + "\t" + affiliatedFactor + "\t" + str(cooCountHash[primaryFactor][affiliatedFactor]) + "\n" )
                else:
                    fo.write( primaryFactor + "\t" + affiliatedFactor + "\t" + "0" + "\n" )
    # done
    #
#
def tab2matrix(tabFile, matFile):
    # table  to matrix
    md = pandas.read_table(tabFile,sep="\t")
    md = md.pivot("primaryFactor","affiliatedFactor","percentage")
    md.to_csv(matFile, sep="\t")
#


inputDIR    = "example_data/expressedNrTfbsInSummit201bp/asTable/"
tabCountDIR = "example_data/observedOccurrence/countAsTable/"
tabFractDIR = "example_data/observedOccurrence/fractionAsTable/"
matFractDIR = "example_data/observedOccurrence/fractionAsMatrix/"

sampleNameList = []

os.system("mkdir -p " + tabCountDIR)
os.system("mkdir -p " + tabFractDIR)
os.system("mkdir -p " + matFractDIR)


filelist = os.listdir(inputDIR)
infoLine("Prepare jobs to generate cooccurance data")
mpPool = multiprocessing.Pool()
for myfile in tqdm(filelist):
    infile = inputDIR + "/" + myfile
    row = myfile.split(".")
    sampleName = row[0]
    sampleNameList.append(sampleName)
    
    tabCountFile= tabCountDIR+ "/" + sampleName + ".pairedColocalizedSummit.tab"
    tabFractFile= tabFractDIR+ "/" + sampleName + ".pairedColocalizedSummit.tab"
    #mpPool.apply_async( processData, args=(infile, tabCountFile, tabFractFile) )
    #
#
infoLine("close and wait pool")
mpPool.close()
mpPool.join()


infoLine("Convert table to matrix")
mpPool = multiprocessing.Pool()
for sampleName in sampleNameList:
    tabFile = tabFractDIR + "/" + sampleName + ".pairedColocalizedSummit.tab"
    matFile = matFractDIR + "/" + sampleName + ".pairedColocalizedSummit.mat"
    mpPool.apply_async( tab2matrix, args=(tabFile, matFile) )
#
infoLine("close and wait pool")
mpPool.close()
mpPool.join()


infoLine("Done!")



