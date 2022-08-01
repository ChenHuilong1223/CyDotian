import os, time, argparse, codecs, sys
sys.stdout = codecs.getwriter("utf-8")(sys.stdout.detach())

parser = argparse.ArgumentParser(description='Batch calculation of internal repetition densities')
parser.add_argument('-l', '--inputRepeatLen', type=int, metavar='length', required=True, help='Input repeat length threshold')
parser.add_argument('-r', '--inputResultFolderPath', metavar='result', required=True, help='Input result folder path')
args = parser.parse_args()

startTime = time.time()
state = 0
try:
    repeatLen = args.inputRepeatLen
    importFolder = args.inputResultFolderPath
    importFolderPath = os.path.abspath(importFolder)

    def mkdir(path):
        folder = os.path.exists(path)
        if not folder:
            os.makedirs(path)
        else:
            pass

    def calculateDirectRepeatDensity(directFilePath, repeatLen, length):
        directPosFile = open(directFilePath, 'r', encoding='utf-8')
        directCount = 0
        for line in directPosFile:
            lineList = line.strip('\n').split('\t')
            if int(lineList[4]) >= repeatLen:
                directCount += 1
        directDensity = (directCount - 1) / (length * length - length)
        directDensity = '{:6f}'.format(directDensity)
        directPosFile.close()
        return directDensity
    def calculateInvertedRepeatDensity(invertedFilePath, repeatLen, length):
        invertedPosFile = open(invertedFilePath, 'r', encoding='utf-8')
        invertedCount = 0
        for line in invertedPosFile:
            lineList = line.strip('\n').split('\t')
            if int(lineList[4]) >= repeatLen:
                invertedCount += 1
        invertedDensity = invertedCount / (length * length)
        invertedDensity = '{:6f}'.format(invertedDensity)
        invertedPosFile.close()
        return invertedDensity
    def calculateRepeatDensity(directFilePath, invertedFilePath, repeatLen, length):
        directPosFile = open(directFilePath, 'r', encoding='utf-8')
        directCount = 0
        for line in directPosFile:
            lineList = line.strip('\n').split('\t')
            if int(lineList[4]) >= repeatLen:
                directCount += 1

        invertedPosFile = open(invertedFilePath, 'r', encoding='utf-8')
        invertedCount = 0
        for line in invertedPosFile:
            lineList = line.strip('\n').split('\t')
            if int(lineList[4]) >= repeatLen:
                invertedCount += 1

        directInvertedDensity = (directCount - 1 + invertedCount) / (length * length)
        directInvertedDensity = '{:6f}'.format(directInvertedDensity)

        directPosFile.close()
        invertedPosFile.close()
        return directInvertedDensity

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

        exportFolder = importFolderPath + '/' + folderName + '/positions_original_densities'
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

        compositionDic = {}
        composition = ''
        positionsOriginalFileGeneNameDic = {}
        for fileName in positionsOriginalFileList:
            name = fileName.split('_positions_')[0]
            matchStr = fileName.split('_positions_')[1]
            compositionDic[matchStr] = ''
            positionsOriginalFileGeneNameDic[name] = ''
            
        positionsOriginalFileGeneNameList = list(positionsOriginalFileGeneNameDic.keys())

        if 'direct.txt' in list(compositionDic.keys()):
            composition += '0'
        if 'inverted.txt' in list(compositionDic.keys()):
            composition += '1'

        directDensityFile, invertedDensityFile, densityFile = None, None, None
        if '0' == composition:
            directDensityFile = open(exportFolderPath + '/' + folderName + '_direct_densities.txt', 'w', encoding='utf-8')
        if '1' == composition:
            invertedDensityFile = open(exportFolderPath + '/' + folderName + '_inverted_densities.txt', 'w', encoding='utf-8')
        if '01' == composition:
            densityFile = open(exportFolderPath + '/' + folderName + '_total_densities.txt', 'w', encoding='utf-8')
            directDensityFile = open(exportFolderPath + '/' + folderName + '_direct_densities.txt', 'w', encoding='utf-8')
            invertedDensityFile = open(exportFolderPath + '/' + folderName + '_inverted_densities.txt', 'w', encoding='utf-8')

        for name in positionsOriginalFileGeneNameList:
            try:
                length = int(seqLengthDic[name])

                if '0' == composition:
                    directFileName = name + '_positions_direct.txt'
                    if directFileName in positionsOriginalFileList:
                        directFilePath = os.path.join(positionsOriginalFolderPath, directFileName)
                        directRepeatDensity = calculateDirectRepeatDensity(directFilePath, repeatLen, length)
                        directDensityFile.write(name+'\t'+ str(directRepeatDensity) + '\n')
                elif '1' == composition:
                    invertedFileName = name + '_positions_inverted.txt'
                    if invertedFileName in positionsOriginalFileList:
                        invertedFilePath = os.path.join(positionsOriginalFolderPath, invertedFileName)
                        invertedRepeatDensity = calculateInvertedRepeatDensity(invertedFilePath, repeatLen, length)
                        invertedDensityFile.write(name + '\t' + str(invertedRepeatDensity) + '\n')
                elif '01' == composition:
                    directFileName = name + '_positions_direct.txt'
                    directFilePath = os.path.join(positionsOriginalFolderPath, directFileName)
                    directRepeatDensity = calculateDirectRepeatDensity(directFilePath, repeatLen, length)
                    directDensityFile.write(name + '\t' + str(directRepeatDensity) + '\n')

                    invertedFileName = name + '_positions_inverted.txt'
                    invertedFilePath = os.path.join(positionsOriginalFolderPath, invertedFileName)
                    invertedRepeatDensity = calculateInvertedRepeatDensity(invertedFilePath, repeatLen, length)
                    invertedDensityFile.write(name + '\t' + str(invertedRepeatDensity) + '\n')

                    totalRepeatDesity = calculateRepeatDensity(directFilePath, invertedFilePath, repeatLen, length)
                    densityFile.write(name + '\t' + str(totalRepeatDesity) + '\n')
            except BaseException as e:
                print(name, e, e.__traceback__.tb_lineno, sep='***')
                failFile.write(name + '\n')

        logFile.close()
        failFile.close()
        if '0' == composition:
            directDensityFile.close()
        if '1' == composition:
            invertedDensityFile.close()
        if '01' == composition:
            densityFile.close()
            directDensityFile.close()
            invertedDensityFile.close()
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

# created by Huilong Chen, May 24, 2022!
# revised by Huilong Chen, May 27, 2022!
# revised by Huilong Chen, July 26, 2022! Optimize.
# revised by Huilong Chen, August 1, 2022! Optimize.
