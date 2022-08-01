import os, re, time, argparse, codecs, sys
import matplotlib.pyplot as plt
sys.stdout = codecs.getwriter("utf-8")(sys.stdout.detach())

parser = argparse.ArgumentParser(description='Batch draw dotplots')
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
    def directDotplot(filePath, repeatLen):
        directPosFile = open(filePath, 'r', encoding='utf-8')
        x, y = [], []
        for line in directPosFile:
            lineList = line.strip('\n').split('\t')
            if int(lineList[4]) >= repeatLen:
                for i in range(int(lineList[0]), int(lineList[0]) + int(lineList[4])):
                    y.append(i)
                for j in range(int(lineList[2]), int(lineList[2]) + int(lineList[4])):
                    x.append(j)
        directPosFile.close()
        return x, y
    def invertedDotplot(filePath, repeatLen):
        invertedPosFile = open(filePath, 'r', encoding='utf-8')
        x, y = [], []
        for line in invertedPosFile:
            lineList = line.strip('\n').split('\t')
            if int(lineList[4]) >= repeatLen:
                for i in range(int(lineList[0]), int(lineList[0]) + int(lineList[4])):
                    y.append(i)
                for j in range(int(lineList[2]), int(lineList[2]) - int(lineList[4]), -1):
                    x.append(j)
        invertedPosFile.close()
        return x, y
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

        exportFolder = importFolderPath+ '/' + folderName + '/positions_original_dotplots'
        mkdir(exportFolder)
        exportFolderPath = os.path.abspath(exportFolder)

        logFile = open(exportFolderPath + '/log.txt', 'w', encoding='utf-8')
        logFile.write('# RepeatLen >= ' + str(repeatLen) + '\n')
        logFile.write('# The names of files that failed to be analysed are recorded in file failure_files.txt.\n')
        logFile.write('Time: ' + time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))

        failFile = open(exportFolderPath + '/failure_files.txt', 'w', encoding='utf-8')
        seqLengthFilePath = os.path.join(positionsOriginalFolderPath, 'all_sequences_length.txt')
        seqLengthDic = creatSeqLengthDic(seqLengthFilePath)
        while 'all_sequences_length.txt' in positionsOriginalFileList:
            positionsOriginalFileList.remove('all_sequences_length.txt')

        for fileName in positionsOriginalFileList:
            try:
                filePath = os.path.join(positionsOriginalFolderPath, fileName)
                if fileName.endswith('_positions_direct.txt'):
                    name = re.findall("(.*?)_positions_direct.txt", fileName)[0]
                    length = int(seqLengthDic[name])
                    fig, ax = plt.subplots()
                    ax.xaxis.set_ticks_position('top')
                    ax.invert_yaxis()
                    ax.set_aspect(1)
                    x, y = directDotplot(filePath,repeatLen)
                    ax.scatter(x, y, marker='o', c='k', s=0.5, linewidths=0, edgecolors=None, alpha=1)
                    ax.set_xlim(0, length+1)
                    ax.set_ylim(length+1, 0)
                    ax.set_xlabel(name, fontdict=fontDict)
                    ax.xaxis.set_label_position('top')
                    ax.set_ylabel(name, fontdict=fontDict)
                    fig.savefig(exportFolderPath + '/' + name + '_direct.png', format='png', dpi=500)
                    fig.savefig(exportFolderPath + '/' + name + '_direct.pdf', format='pdf')
                    plt.clf()
                elif fileName.endswith('_positions_inverted.txt'):
                    name = re.findall("(.*?)_positions_inverted.txt", fileName)[0]
                    length = int(seqLengthDic[name])
                    fig, ax = plt.subplots()
                    ax.xaxis.set_ticks_position('top')
                    ax.invert_yaxis()
                    ax.set_aspect(1)
                    x, y = invertedDotplot(filePath,repeatLen)
                    ax.scatter(x, y, marker='o', c='k', s=0.5, linewidths=0, edgecolors=None, alpha=1)
                    ax.set_xlim(0, length+1)
                    ax.set_ylim(length+1, 0)
                    ax.set_xlabel(name, fontdict=fontDict)
                    ax.xaxis.set_label_position('top')
                    ax.set_ylabel(name, fontdict=fontDict)
                    fig.savefig(exportFolderPath + '/' + name + '_inverted.png', format='png', dpi=500)
                    fig.savefig(exportFolderPath + '/' + name + '_inverted.pdf', format='pdf')
                    plt.clf()
                plt.close("all")
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

# created by Huilong Chen, May 22, 2022!
# revised by Huilong Chen, May 27, 2022!
# revised by Huilong Chen, July 26, 2022! Optimize.
# revised by Huilong Chen, August 1, 2022! Optimize.
