import os, re, time, argparse, codecs, sys

sys.stdout = codecs.getwriter("utf-8")(sys.stdout.detach())

parser = argparse.ArgumentParser(description='Batch extract repeat sequences')
parser.add_argument('-l', '--inputRepeatLen', type=int, metavar='length', required=True, help='Input repeat length threshold')
parser.add_argument('-r', '--inputResultFolderPath', metavar='result', required=True, help='Input result folder path')
parser.add_argument('-s', '--inputSequenceFolderPath', metavar='sequence', required=True, help='Input sequence folder path')
args = parser.parse_args()

startTime = time.time()
state = 0
try:
    repeatLen = args.inputRepeatLen
    importFolder = args.inputResultFolderPath
    importFolderPath = os.path.abspath(importFolder)

    importSequenceFolder = args.inputSequenceFolderPath
    importSequenceFileList = os.listdir(importSequenceFolder)
    importSequenceFolderPath = os.path.abspath(importSequenceFolder)


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

    # 这里可以判断importSequenceFileList和folderList是否有差集，否则抛出异常!
    if [] != list(set(folderList) - set(importSequenceFileList)):
        print('The corresponding fasta file is missing in the input folder containing fasta files, please check!')
        sys.exit(-1)


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


    for folderName in folderList:
        sequencefilePath = os.path.join(importSequenceFolderPath, folderName)
        chlFasta = creatFastaDict(sequencefilePath)

        positionsFolder = importFolderPath + '/' + folderName + '/positions'
        positionsFileList = os.listdir(positionsFolder)
        positionsFolderPath = os.path.abspath(positionsFolder)

        exportFolder = importFolderPath + '/' + folderName + '/positions_repeat_sequences'
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
                    directRepeatSequenceFilePath = exportFolderPath + '/' + name + '_direct_repeat_sequences.txt'
                    with open(filePath, 'r', encoding='utf-8') as positionFile:
                        for line in positionFile:
                            lineList = line.strip('\n').split('\t')
                            if int(lineList[4]) >= repeatLen:
                                directRepeatSequenceFile = open(directRepeatSequenceFilePath, 'w', encoding='utf-8') # 有大于指定长度的结果才创建吧，不然用户会说无法判断是没有，还是程序输出错误。
                            break
                    # posTable = pd.read_csv(filePath, encoding='utf-8', sep='\t', header=None)
                    # for chl in posTable.values:
                    # positionFile = open(filePath, 'r', encoding='utf-8')
                    with open(filePath, 'r', encoding='utf-8') as positionFile:
                        for line in positionFile:
                            lineList = line.strip('\n').split('\t')
                            if int(lineList[4]) >= repeatLen:
                                print(line, end='', file=directRepeatSequenceFile)
                                print('vertical', file=directRepeatSequenceFile)
                                print(chlFasta[name][int(lineList[0]) - 1:int(lineList[1])], file=directRepeatSequenceFile)
                                print(chlFasta[name][int(lineList[2]) - 1:int(lineList[3])], file=directRepeatSequenceFile)
                                print('horizontal\n', file=directRepeatSequenceFile)

                    with open(filePath, 'r', encoding='utf-8') as positionFile:
                        for line in positionFile:
                            lineList = line.strip('\n').split('\t')
                            if int(lineList[4]) >= repeatLen:
                                directRepeatSequenceFile.close()
                            break
                elif fileName.endswith('_positions_inverted.txt'):
                    name = re.findall("(.*?)_positions_inverted.txt", fileName)[0]
                    invertedRepeatSequenceFilePath = exportFolderPath + '/' + name + '_inverted_repeat_sequences.txt'
                    with open(filePath, 'r', encoding='utf-8') as positionFile:
                        for line in positionFile:
                            lineList = line.strip('\n').split('\t')
                            if int(lineList[4]) >= repeatLen:
                                invertedRepeatSequenceFile = open(invertedRepeatSequenceFilePath, 'w', encoding='utf-8') # 有大于指定长度的结果才创建吧，不然用户会说无法判断是没有，还是程序输出错误。
                            break

                    with open(filePath, 'r', encoding='utf-8') as positionFile:
                        for line in positionFile:
                            lineList = line.strip('\n').split('\t')
                            if int(lineList[4]) >= repeatLen:
                                print(line, end='', file=invertedRepeatSequenceFile)
                                print('vertical', file=invertedRepeatSequenceFile)
                                print(chlFasta[name][int(lineList[0]) - 1:int(lineList[1])], file=invertedRepeatSequenceFile)
                                print(chlFasta[name][int(lineList[3]) - 1:int(lineList[2])], file=invertedRepeatSequenceFile)
                                print('horizontal\n', file=invertedRepeatSequenceFile)

                    with open(filePath, 'r', encoding='utf-8') as positionFile:
                        for line in positionFile:
                            lineList = line.strip('\n').split('\t')
                            if int(lineList[4]) >= repeatLen:
                                invertedRepeatSequenceFile.close()
                            break
                elif fileName.endswith('_positions_reverse_complement.txt'):
                    name = re.findall("(.*?)_positions_reverse_complement.txt", fileName)[0]
                    reverseComplementRepeatSequenceFilePath = exportFolderPath + '/' + name + '_reverse_complement_sequences.txt'
                    with open(filePath, 'r', encoding='utf-8') as positionFile:
                        for line in positionFile:
                            lineList = line.strip('\n').split('\t')
                            if int(lineList[4]) >= repeatLen:
                                reverseComplementRepeatSequenceFile = open(reverseComplementRepeatSequenceFilePath, 'w',
                                                                  encoding='utf-8')  # 有大于指定长度的结果才创建吧，不然用户会说无法判断是没有，还是程序输出错误。
                            break

                    with open(filePath, 'r', encoding='utf-8') as positionFile:
                        for line in positionFile:
                            lineList = line.strip('\n').split('\t')
                            if int(lineList[4]) >= repeatLen:
                                print(line, end='', file=reverseComplementRepeatSequenceFile)
                                print('vertical', file=reverseComplementRepeatSequenceFile)
                                print(chlFasta[name][int(lineList[0]) - 1:int(lineList[1])],
                                      file=reverseComplementRepeatSequenceFile)
                                print(chlFasta[name][int(lineList[3]) - 1:int(lineList[2])],
                                      file=reverseComplementRepeatSequenceFile)
                                print('horizontal\n', file=reverseComplementRepeatSequenceFile)

                    with open(filePath, 'r', encoding='utf-8') as positionFile:
                        for line in positionFile:
                            lineList = line.strip('\n').split('\t')
                            if int(lineList[4]) >= repeatLen:
                                reverseComplementRepeatSequenceFile.close()
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
# revised by Huilong Chen, August 1, 2022! Optimize.
