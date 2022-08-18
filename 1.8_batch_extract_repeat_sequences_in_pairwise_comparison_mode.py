import os, re, time, argparse, codecs, sys

sys.stdout = codecs.getwriter("utf-8")(sys.stdout.detach())

parser = argparse.ArgumentParser(description='Batch extract repeat sequences')
parser.add_argument('-l', '--inputRepeatLen', type=int, metavar='length', required=True, help='Input repeat length threshold')
parser.add_argument('-r', '--inputResultFolderPath', metavar='result', required=True, help='Input result folder path')
parser.add_argument('-s1', '--inputVerticalSequenceFilePath', metavar='vertical', required=True, help='Input vertical sequence file path')
parser.add_argument('-s2', '--inputHorizontalSequenceFilePath', metavar='horizontal', required=True, help='Input horizontal sequence file path')
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

    folderList = []
    for root, dirs, files in os.walk(importFolder):
        folderList = dirs
        break

    def creatFastaDict(sequencefilePath):
        with open(sequencefilePath, 'r', encoding='utf-8') as fastaFile:
            chlFasta = {}
            gene = seq = ''
            for row in fastaFile:
                row = row.strip('\n')
                if row.startswith('>'):
                    if gene != '' and seq != '':
                        chlFasta[gene] = seq.upper()
                    gene = row.replace('>', '')
                    seq = ''
                else:
                    seq += row
            chlFasta[gene] = seq.upper()
        return chlFasta

    chlFasta1 = creatFastaDict(args.inputVerticalSequenceFilePath)
    chlFasta2 = creatFastaDict(args.inputHorizontalSequenceFilePath)

    positionsFolder = importFolderPath + '/positions_original'
    positionsFileList = os.listdir(positionsFolder)
    positionsFolderPath = os.path.abspath(positionsFolder)

    exportFolder = importFolderPath + '/positions_original_similar_sequences'
    mkdir(exportFolder)
    exportFolderPath = os.path.abspath(exportFolder)

    logFile = open(exportFolderPath + '/log.txt', 'w', encoding='utf-8')
    logFile.write('# RepeatLen >= ' + str(repeatLen) + '\n')
    logFile.write('# The names of files that failed to be analysed are recorded in file failure_files.txt.\n')
    logFile.write('Time: ' + time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))

    failFile = open(exportFolderPath + '/failure_files.txt', 'w', encoding='utf-8')

    for fileName in positionsFileList:
        try:
            filePath = os.path.join(positionsFolderPath, fileName)
            if fileName.endswith('_positions_direct.txt'):
                name = re.findall("(.*?)_positions_direct.txt", fileName)[0]
                name1 = name.split('_VS_')[0]
                name2 = name.split('_VS_')[1]
                directSimilarSequenceFilePath = exportFolderPath + '/' + name + '_direct_similar_sequences.txt'
                with open(filePath, 'r', encoding='utf-8') as positionFile:
                    for line in positionFile:
                        lineList = line.strip('\n').split('\t')
                        if int(lineList[4]) >= repeatLen:
                            directSimilarSequenceFile = open(directSimilarSequenceFilePath, 'w', encoding='utf-8')
                        break

                if name1 != name2:
                    with open(filePath, 'r', encoding='utf-8') as positionFile:
                        for line in positionFile:
                            lineList = line.strip('\n').split('\t')
                            if int(lineList[4]) >= repeatLen:
                                print(line, end='', file=directSimilarSequenceFile)
                                print(name1, 'vertical', sep=', ', file=directSimilarSequenceFile)
                                print(chlFasta1[name1][int(lineList[0]) - 1:int(lineList[1])], file=directSimilarSequenceFile)
                                print(chlFasta2[name2][int(lineList[2]) - 1:int(lineList[3])], file=directSimilarSequenceFile)
                                print(name2, 'horizontal\n', sep=', ', file=directSimilarSequenceFile)
                else:
                    with open(filePath, 'r', encoding='utf-8') as positionFile:
                        for line in positionFile:
                            lineList = line.strip('\n').split('\t')
                            if int(lineList[4]) >= repeatLen and int(lineList[2]) > int(lineList[0]):
                                print(line, end='', file=directSimilarSequenceFile)
                                print(name1, 'vertical', sep=', ', file=directSimilarSequenceFile)
                                print(chlFasta1[name1][int(lineList[0]) - 1:int(lineList[1])], file=directSimilarSequenceFile)
                                print(chlFasta2[name2][int(lineList[2]) - 1:int(lineList[3])], file=directSimilarSequenceFile)
                                print(name2, 'horizontal\n', sep=', ', file=directSimilarSequenceFile)

                with open(filePath, 'r', encoding='utf-8') as positionFile:
                    for line in positionFile:
                        lineList = line.strip('\n').split('\t')
                        if int(lineList[4]) >= repeatLen:
                            directSimilarSequenceFile.close()
                        break
            elif fileName.endswith('_positions_inverted.txt'):
                name = re.findall("(.*?)_positions_inverted.txt", fileName)[0]
                name1 = name.split('_VS_')[0]
                name2 = name.split('_VS_')[1]
                invertedSimilarSequenceFilePath = exportFolderPath + '/' + name + '_inverted_similar_sequences.txt'
                with open(filePath, 'r', encoding='utf-8') as positionFile:
                    for line in positionFile:
                        lineList = line.strip('\n').split('\t')
                        if int(lineList[4]) >= repeatLen:
                            invertedSimilarSequenceFile = open(invertedSimilarSequenceFilePath, 'w', encoding='utf-8')
                        break

                if name1 != name2:
                    with open(filePath, 'r', encoding='utf-8') as positionFile:
                        for line in positionFile:
                            lineList = line.strip('\n').split('\t')
                            if int(lineList[4]) >= repeatLen:
                                print(line, end='', file=invertedSimilarSequenceFile)
                                print(name1, 'vertical', sep=', ', file=invertedSimilarSequenceFile)
                                print(chlFasta1[name1][int(lineList[0]) - 1:int(lineList[1])], file=invertedSimilarSequenceFile)
                                print(chlFasta2[name2][int(lineList[3]) - 1:int(lineList[2])], file=invertedSimilarSequenceFile)
                                print(name2, 'horizontal\n', sep=', ', file=invertedSimilarSequenceFile)
                else:
                    with open(filePath, 'r', encoding='utf-8') as positionFile:
                        for line in positionFile:
                            lineList = line.strip('\n').split('\t')
                            if int(lineList[4]) >= repeatLen and int(lineList[2]) > int(lineList[0]):
                                print(line, end='', file=invertedSimilarSequenceFile)
                                print(name1, 'vertical', sep=', ', file=invertedSimilarSequenceFile)
                                print(chlFasta1[name1][int(lineList[0]) - 1:int(lineList[1])], file=invertedSimilarSequenceFile)
                                print(chlFasta2[name2][int(lineList[3]) - 1:int(lineList[2])], file=invertedSimilarSequenceFile)
                                print(name2, 'horizontal\n', sep=', ', file=invertedSimilarSequenceFile)

                with open(filePath, 'r', encoding='utf-8') as positionFile:
                    for line in positionFile:
                        lineList = line.strip('\n').split('\t')
                        if int(lineList[4]) >= repeatLen:
                            invertedSimilarSequenceFile.close()
                        break
        except BaseException as e:
            print(fileName, e, e.__traceback__.tb_lineno, sep='***')
            failFile.write(fileName + '\n')

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

# created by Huilong Chen, July 28, 2022!
