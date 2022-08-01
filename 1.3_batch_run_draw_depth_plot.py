import os, re, time, argparse, codecs, sys
import matplotlib.pyplot as plt
import numpy as np
sys.stdout = codecs.getwriter("utf-8")(sys.stdout.detach())

parser = argparse.ArgumentParser(description='Batch draw depth plots')
parser.add_argument('-l', '--inputRepeatLen', type=int, metavar='length', required=True, help='Input repeat length threshold')
parser.add_argument('-r', '--inputResultFolderPath', metavar='result', required=True, help='Input result folder path')
args = parser.parse_args()

startTime = time.time()
state = 0
try:
    plt.rcParams['font.family'] = ["Times New Roman"]
    fontDict = {"size": 13, "color": "k", 'family': 'Times New Roman'}
    repeatLen = args.inputRepeatLen
    importFolder = args.inputResultFolderPath
    importFolderPath = os.path.abspath(importFolder)

    def mkdir(path):
        folder = os.path.exists(path)
        if not folder:
            os.makedirs(path)
        else:
            pass
    def directDepth(directFilePath, repeatLen):
        directPosFile = open(directFilePath, 'r', encoding='utf-8')
        directDepthArray = np.full(int(length), -1, dtype=int, order='C')
        for line in directPosFile:
            lineList = line.strip('\n').split('\t')
            if int(lineList[4]) >= repeatLen:
                for j in range(int(lineList[2]), int(lineList[2]) + int(lineList[4])):
                    directDepthArray[j - 1] += 1
        directPosFile.close()
        return directDepthArray
    def invertedDepth(invertedFilePath, repeatLen):
        invertedPosFile = open(invertedFilePath, 'r', encoding='utf-8')
        invertedDepthArray = np.zeros(int(length), dtype=int, order='C')
        for line in invertedPosFile:
            lineList = line.strip('\n').split('\t')
            if int(lineList[4]) >= repeatLen:
                for j in range(int(lineList[2]), int(lineList[2]) - int(lineList[4]), -1):
                    invertedDepthArray[j - 1] -= 1
        invertedPosFile.close()
        return invertedDepthArray
    def creatSeqLengthDic(seqLengthFilePath):
        seqLengthFile = open(seqLengthFilePath, 'r', encoding='utf-8')
        seqLengthDic = {}
        for line in seqLengthFile:
            lineList = line.strip('\n').rsplit('\t',1)
            seqLengthDic[lineList[0]] = lineList[1]
        seqLengthFile.close()
        return seqLengthDic

    folderList = []
    for root, dirs, files in os.walk(importFolder):
        folderList = dirs
        break
    for folderName in folderList:
        positionsOriginalFolder = importFolderPath + '/' + folderName + '/positions_original'
        positionsOriginalFileList = os.listdir(positionsOriginalFolder)
        positionsOriginalFolderPath = os.path.abspath(positionsOriginalFolder)

        exportFolder = importFolderPath+ '/' + folderName + '/positions_original_depths'
        mkdir(exportFolder)
        exportFolderPath = os.path.abspath(exportFolder)

        logFile = open(exportFolderPath + '/log.txt', 'w', encoding='utf-8')
        logFile.write('# RepeatLen >= ' + str(repeatLen) + '\n')
        logFile.write('# The gene names of files that failed to be analysed are recorded in file failure_files.txt.\n')
        logFile.write('Time: ' + time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))

        failFile = open(exportFolderPath + '/failure_files.txt', 'w', encoding='utf-8')

        seqLengthFilePath = os.path.join(positionsOriginalFolderPath, 'all_sequences_length.txt')
        seqLengthDic = creatSeqLengthDic(seqLengthFilePath)
        while 'all_sequences_length.txt' in positionsOriginalFileList:
            positionsOriginalFileList.remove('all_sequences_length.txt')

        positionsOriginalFileGeneNameDic = {}
        for fileName in positionsOriginalFileList:
            name = fileName.split('_positions_')[0]
            positionsOriginalFileGeneNameDic[name] = ''
        positionsOriginalFileGeneNameList = list(positionsOriginalFileGeneNameDic.keys())
        for name in positionsOriginalFileGeneNameList:
            try:
                length = int(seqLengthDic[name])

                fig, ax = plt.subplots()
                ax.set_xlim(0, length + 1)
                ax.set_xlabel('Position', fontdict=fontDict)
                ax.set_ylabel('Depth', fontdict=fontDict)
                ax.set_title(name, fontdict=fontDict)
                x = [j for j in range(1, length + 1)]

                directFileName = name + '_positions_direct.txt'
                if directFileName in positionsOriginalFileList:
                    directFilePath = os.path.join(positionsOriginalFolderPath, directFileName)
                    directDepthArray = directDepth(directFilePath,repeatLen)
                    ax.bar(x, height=directDepthArray, width=1.0, alpha=0.7, color='blue')

                    depthFile = open(exportFolderPath + '/' + name + '_direct_depths.txt', 'w')
                    for a in directDepthArray:
                        depthFile.write(str(a) + '\n')
                    depthFile.close()

                invertedFileName = name + '_positions_inverted.txt'
                if invertedFileName in positionsOriginalFileList:
                    invertedFilePath = os.path.join(positionsOriginalFolderPath, invertedFileName)
                    invertedDepthArray = invertedDepth(invertedFilePath, repeatLen)
                    ax.bar(x, height=invertedDepthArray, width=1.0, alpha=0.7, color='red')

                    depthFile = open(exportFolderPath + '/' + name + '_inverted_depths.txt', 'w')
                    for a in invertedDepthArray:
                        depthFile.write(str(abs(a)) + '\n')
                    depthFile.close()

                fig.savefig(exportFolderPath + '/' + name + '_depth.png', format='png', dpi=500)
                fig.savefig(exportFolderPath + '/' + name + '_depth.pdf', format='pdf')
                plt.clf()
                plt.close("all")
            except BaseException as e:
                print(name, e, e.__traceback__.tb_lineno, sep='***')
                failFile.write(name + '\n')

        logFile.close()
        failFile.close()
except BaseException as e:
    state = -1
    print(e, e.__traceback__.tb_lineno)

if 0 == state:
    print('Congratulations, the script worked and finished successfully!')
elif -1 == state:
    print('Sadly, the script did not complete properly, please check the output log, resolve the problem, and try again!')

endTime = time.time()
runTime = round(endTime - startTime)
hour = runTime//3600
minute = (runTime-3600*hour)//60
second = runTime-3600*hour-60*minute
print(f'The program running time: {hour}hour(s) {minute}minute(s) {second}second(s).')

# created by Huilong Chen, May 23, 2022!
# revised by Huilong Chen, May 27, 2022!
# revised by Huilong Chen, July 26, 2022! Optimize.
# revised by Huilong Chen, August 1, 2022! Optimize.
