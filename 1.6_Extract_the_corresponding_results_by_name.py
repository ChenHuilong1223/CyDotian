import os, re, time, argparse, codecs, sys, shutil

sys.stdout = codecs.getwriter("utf-8")(sys.stdout.detach())

parser = argparse.ArgumentParser(description='Batch extract repeat sequences')
parser.add_argument('-n', '--inputNameFilePath', metavar='name', required=True, help='Input name file path')
parser.add_argument('-r', '--inputOneResultFolderPath', metavar='result', required=True, help='Input one result folder path')
parser.add_argument('-o', '--outputOneFolderPath', metavar='output', required=True, help='Output one folder path')
args = parser.parse_args()

startTime = time.time()
state = 0
try:
    inputNameFilePath = os.path.abspath(args.inputNameFilePath)
    nameFileName = os.path.basename(inputNameFilePath)
    with open(inputNameFilePath, 'r', encoding='utf-8') as nameFile:
        nameList = []
        for line in nameFile:
            name = line.strip('\n')
            nameList.append(name)

    importFolder = args.inputOneResultFolderPath
    importFolderPath = os.path.abspath(importFolder)

    exportFolder = args.outputOneFolderPath
    def mkdir(path):
        folder = os.path.exists(path)
        if not folder:
            os.makedirs(path)
        else:
            pass

    mkdir(exportFolder)
    exportFolderPath = os.path.abspath(exportFolder)

    folderList = []
    for root, dirs, files in os.walk(importFolder):
        folderList = dirs
        break

    def copyFile(sourceFilePath, destinationFolderPath):  # 复制函数
        if not os.path.isfile(sourceFilePath):
            print("{} not exist!".format(sourceFilePath))
        else:
            fileFolderPath, fileName = os.path.split(sourceFilePath)  # 分离文件名和路径
            if not os.path.exists(destinationFolderPath):
                os.makedirs(destinationFolderPath)  # 创建路径
            shutil.copy(sourceFilePath, destinationFolderPath + fileName)  # 复制文件

    logFile = open(exportFolderPath + '/log.txt', 'w', encoding='utf-8')
    logFile.write('Time: ' + time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    logFile.close()

    for folderName in folderList:
        resultFolder = importFolderPath + '/' + folderName
        resultFileList = os.listdir(resultFolder)
        if 'positions' == folderName or 'positions_original' == folderName:
            for resultFileName in resultFileList:
                if '_positions_' in resultFileName:
                    name = resultFileName.split('_positions_')[0]
                    if name in nameList:
                        resultFilePath = resultFolder + '/' + resultFileName
                        destinationFolderPath = exportFolderPath + '/' + folderName
                        mkdir(destinationFolderPath)
                        copyFile(resultFilePath, destinationFolderPath + '/')
        elif 'positions_original_dotplots' == folderName:
            for resultFileName in resultFileList:
                # if resultFileName.endswith('_direct.pdf') or resultFileName.endswith('_direct.png') or resultFileName.endswith('_inverted.pdf') or resultFileName.endswith('_inverted.png'):
                if resultFileName.endswith('_direct.pdf'):
                    name = re.findall("(.*?)_direct.pdf", resultFileName)[0]
                    if name in nameList:
                        resultFilePath = resultFolder + '/' + resultFileName
                        destinationFolderPath = exportFolderPath + '/' + folderName
                        mkdir(destinationFolderPath)
                        copyFile(resultFilePath, destinationFolderPath + '/')
                elif resultFileName.endswith('_direct.png'):
                    name = re.findall("(.*?)_direct.png", resultFileName)[0]
                    if name in nameList:
                        resultFilePath = resultFolder + '/' + resultFileName
                        destinationFolderPath = exportFolderPath + '/' + folderName
                        mkdir(destinationFolderPath)
                        copyFile(resultFilePath, destinationFolderPath + '/')
                elif resultFileName.endswith('_inverted.pdf'):
                    name = re.findall("(.*?)_inverted.pdf", resultFileName)[0]
                    if name in nameList:
                        resultFilePath = resultFolder + '/' + resultFileName
                        destinationFolderPath = exportFolderPath + '/' + folderName
                        mkdir(destinationFolderPath)
                        copyFile(resultFilePath, destinationFolderPath + '/')
                elif resultFileName.endswith('_inverted.png'):
                    name = re.findall("(.*?)_inverted.png", resultFileName)[0]
                    if name in nameList:
                        resultFilePath = resultFolder + '/' + resultFileName
                        destinationFolderPath = exportFolderPath + '/' + folderName
                        mkdir(destinationFolderPath)
                        copyFile(resultFilePath, destinationFolderPath + '/')
        elif 'positions_original_depths' == folderName:
            for resultFileName in resultFileList:
                if resultFileName.endswith('_depth.pdf'):
                    name = re.findall("(.*?)_depth.pdf", resultFileName)[0]
                    if name in nameList:
                        resultFilePath = resultFolder + '/' + resultFileName
                        destinationFolderPath = exportFolderPath + '/' + folderName
                        mkdir(destinationFolderPath)
                        copyFile(resultFilePath, destinationFolderPath + '/')
                elif resultFileName.endswith('_depth.png'):
                    name = re.findall("(.*?)_depth.png", resultFileName)[0]
                    if name in nameList:
                        resultFilePath = resultFolder + '/' + resultFileName
                        destinationFolderPath = exportFolderPath + '/' + folderName
                        mkdir(destinationFolderPath)
                        copyFile(resultFilePath, destinationFolderPath + '/')
                elif resultFileName.endswith('_direct_depths.txt'):
                    name = re.findall("(.*?)_direct_depths.txt", resultFileName)[0]
                    if name in nameList:
                        resultFilePath = resultFolder + '/' + resultFileName
                        destinationFolderPath = exportFolderPath + '/' + folderName
                        mkdir(destinationFolderPath)
                        copyFile(resultFilePath, destinationFolderPath + '/')
                elif resultFileName.endswith('_inverted_depths.txt'):
                    name = re.findall("(.*?)_inverted_depths.txt", resultFileName)[0]
                    if name in nameList:
                        resultFilePath = resultFolder + '/' + resultFileName
                        destinationFolderPath = exportFolderPath + '/' + folderName
                        mkdir(destinationFolderPath)
                        copyFile(resultFilePath, destinationFolderPath + '/')
        elif 'positions_original_density' == folderName:
            for resultFileName in resultFileList:
                if resultFileName.endswith('_direct_densities.txt'):
                    resultFilePath = resultFolder + '/' + resultFileName
                    resultFile = open(resultFilePath, 'r', encoding='utf-8')

                    destinationFilePath = exportFolderPath + '/' + folderName + '/' + nameFileName + '_direct_densities.txt'
                    destinationFile = open(destinationFilePath, 'r', encoding='utf-8')

                    for line in resultFile:
                        lineList = line.strip('\n').split('\t')
                        name = lineList[0]
                        if name in nameList:
                            destinationFile.write(line)

                    resultFile.close()
                    destinationFile.close()
                elif resultFileName.endswith('_inverted_densities.txt'):
                    resultFilePath = resultFolder + '/' + resultFileName
                    resultFile = open(resultFilePath, 'r', encoding='utf-8')

                    destinationFilePath = exportFolderPath + '/' + folderName + '/' + nameFileName + '_inverted_densities.txt'
                    destinationFile = open(destinationFilePath, 'r', encoding='utf-8')

                    for line in resultFile:
                        lineList = line.strip('\n').split('\t')
                        name = lineList[0]
                        if name in nameList:
                            destinationFile.write(line)

                    resultFile.close()
                    destinationFile.close()
                elif resultFileName.endswith('_total_densities.txt'):
                    resultFilePath = resultFolder + '/' + resultFileName
                    resultFile = open(resultFilePath, 'r', encoding='utf-8')

                    destinationFilePath = exportFolderPath + '/' + folderName + '/' + nameFileName + '_total_densities.txt'
                    destinationFile = open(destinationFilePath, 'r', encoding='utf-8')

                    for line in resultFile:
                        lineList = line.strip('\n').split('\t')
                        name = lineList[0]
                        if name in nameList:
                            destinationFile.write(line)

                    resultFile.close()
                    destinationFile.close()
        elif 'positions_repeat_sequences' == folderName:
            for resultFileName in resultFileList:
                if resultFileName.endswith('_direct_repeat_sequences.txt'):
                    name = re.findall("(.*?)_direct_repeat_sequences.txt", resultFileName)[0]
                    if name in nameList:
                        resultFilePath = resultFolder + '/' + resultFileName
                        destinationFolderPath = exportFolderPath + '/' + folderName
                        mkdir(destinationFolderPath)
                        copyFile(resultFilePath, destinationFolderPath + '/')
                elif resultFileName.endswith('_inverted_repeat_sequences.txt'):
                    name = re.findall("(.*?)_inverted_repeat_sequences.txt", resultFileName)[0]
                    if name in nameList:
                        resultFilePath = resultFolder + '/' + resultFileName
                        destinationFolderPath = exportFolderPath + '/' + folderName
                        mkdir(destinationFolderPath)
                        copyFile(resultFilePath, destinationFolderPath + '/')
                elif resultFileName.endswith('_reverse_complement_sequences.txt'):
                    name = re.findall("(.*?)_reverse_complement_sequences.txt", resultFileName)[0]
                    if name in nameList:
                        resultFilePath = resultFolder + '/' + resultFileName
                        destinationFolderPath = exportFolderPath + '/' + folderName
                        mkdir(destinationFolderPath)
                        copyFile(resultFilePath, destinationFolderPath + '/')

except BaseException as e:
    state = -1
    print(e, e.__traceback__.tb_lineno)

if 0 == state:
    print('Congratulations, the script worked and finished successfully!')
elif -1 == state:
    print(
        'Sadly, the script did not complete properly, please check the output log, resolve the problem, and try again!')

endTime = time.time()
runTime = round(endTime - startTime)
hour = runTime//3600
minute = (runTime-3600*hour)//60
second = runTime-3600*hour-60*minute
print(f'The program running time: {hour}hour(s) {minute}minute(s) {second}second(s).')

# created by Huilong Chen, July 30, 2022!
# revised by Huilong Chen, August 1, 2022! Optimize.
