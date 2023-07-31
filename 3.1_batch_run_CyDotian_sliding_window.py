import os, re, time, argparse, codecs, sys, subprocess, signal
import pandas as pd
sys.stdout = codecs.getwriter("utf-8")(sys.stdout.detach())

parser = argparse.ArgumentParser(description="Batch run CyDotian's modules")
parser.add_argument('-c', '--inputConfigFilePath', metavar='config', required=True, help='Input config file path')
parser.add_argument('-s', '--inputSequenceFolderPath', metavar='sequence', required=True, help='Input sequence folder path')
parser.add_argument('-o', '--outputFolderPath', metavar='output', required=True, help='Output folder path')
args = parser.parse_args()

startTime = time.time()
state = 0
try:
    importFolder = args.inputSequenceFolderPath
    importFileList = os.listdir(importFolder)
    importFolderPath = os.path.abspath(importFolder)

    exportFolder = args.outputFolderPath
    def mkdir(path):
        folder = os.path.exists(path)
        if not folder:
            os.makedirs(path)
        else:
            pass
    
    mkdir(exportFolder)
    exportFolderPath = os.path.abspath(exportFolder)

    configFile = open(args.inputConfigFilePath,'r',encoding='utf-8')
    ideSimThr,windowSize,modeList,aminoAcidMatrix,fileType = 0.85, 15, ['0', '1'], 1, 0
    print('=====================')
    for row in configFile:
        if False == row.startswith('#'):
            if '=' in row:
                parameterName = re.findall("(.*?)=", row.strip('\n'))[0]
                matchStr = re.findall("=(.*?)#", row.strip('\n'))[0]
                print(parameterName, matchStr, sep='=')
                parameter = matchStr.replace(' ','')
                if row.startswith('ideSimThr'):
                    ideSimThr = float(parameter)
                elif row.startswith('windowSize'):
                    windowSize = int(parameter)
                elif row.startswith('mode'):
                    modeList = parameter.split(',')
                elif row.startswith('aminoAcidMatrix'):
                    aminoAcidMatrix = int(parameter)
                elif row.startswith('fileType'):
                    if '0' == parameter:
                        fileType = 0
                    elif '1' == parameter:
                        fileType = 1
    configFile.close()
    print('=====================\n')

    def batchExportPosition(chlFasta,exportFolderPositionPath,ideSimThr,windowSize,modeList,aminoAcidMatrix,fileType):
        folder = '/plot_positions'
        customFolderPath = exportFolderPositionPath + folder
        isExists = os.path.exists(customFolderPath)
        if isExists:
            pass
        else:
            os.makedirs(customFolderPath, mode=0o777)
    
        logFile = open(customFolderPath + '/log.txt', 'w',encoding='utf-8')
        if 0 == fileType:
            logFile.write('# DNA' + '\n')
            logFile.write(
                '# parameter: Identity >= ' + str(ideSimThr * 100) + '% ' + 'WindowSize >= ' + str(windowSize) + '\n')
            logFile.write(
                '# header of plot_positions.txt: y  x' + '\n')
        elif 1 == fileType:
            aminoAcidMatrixList = ['BLOSUM45', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90', 'PAM30', 'PAM70', 'PAM250']
            logFile.write('# Amino acid' + '\t' + aminoAcidMatrixList[aminoAcidMatrix] + '\n')
            logFile.write(
                '# parameter: Similarity >= ' + str(ideSimThr * 100) + '% ' + 'WindowSize >= ' + str(
                    windowSize) + '\n')
            logFile.write(
                '# header of positions.txt: y  x' + '\n')
        logFile.write('# The names and lengths of sequences that failed to be analysed are recorded in file failure_sequences.txt.\n')
        logFile.write('Time: ' + time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))

        failFile = open(customFolderPath + '/failure_sequences.txt', 'w',encoding='utf-8')
        failLogFile = open(customFolderPath + '/failure_sequences_error_log.txt', 'w', encoding='utf-8')

        seqLengthFileOriginal = open(customFolderPath + '/all_sequences_length.txt', 'w', encoding='utf-8')

        def removeThreeFiles():
            if os.path.exists('temp.single.input.fasta1.txt'):
                os.remove('temp.single.input.fasta1.txt')
            else:
                print("'temp.single.input.fasta1.txt' does not exist!")
            if os.path.exists('temp.single.input.fasta2.txt'):
                os.remove('temp.single.input.fasta2.txt')
            else:
                print("'temp.single.input.fasta2.txt' does not exist!")
            if os.path.exists('position.txt'):
                os.remove('position.txt')
            else:
                print("'position.txt' does not exist!")

        for name, seq in chlFasta.items():
            try:
                length = len(seq)
                seqLengthFileOriginal.write(name + '\t' + str(length) + '\n')
                if '0' in modeList:
                    tempSingleInputFile1 = open('./temp.single.input.fasta1.txt','w',encoding='utf-8')
                    tempSingleInputFile1.write(seq)
                    tempSingleInputFile1.close()
                    tempSingleInputFile2 = open('./temp.single.input.fasta2.txt','w',encoding='utf-8')
                    tempSingleInputFile2.write(seq)
                    tempSingleInputFile2.close()
                    mode = 0
                    # command = 'cd ./bin/; ./slidingWindow'+' '+str(ideSimThr)+' '+str(windowSize)+' '+str(mode)+' '+str(aminoAcidMatrix)+' '+str(fileType)
                    # os.system(command)
                    command = "./slidingWindow {} {} {} {} {}".format(str(ideSimThr), str(windowSize), str(mode), str(aminoAcidMatrix), str(fileType))
                    process = subprocess.Popen(command, universal_newlines=True, stdin=subprocess.PIPE,
                                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,
                                               preexec_fn=os.setsid)
                    stdoutput, erroutput = process.communicate()  # 这行代码保证调用的子进程结束之后再执行Python脚本中下面的代码。

                    print("mode 0, {}, process.returncode: {}".format(name, process.returncode))
                    if process.returncode:  # 获取进程的返回值。如果进程还没有结束，返回None。根据自我实验，就是等价于Popen.poll()
                        print(name, str(length), 'errOutput: ', erroutput.strip('\n'), sep=', ', end='***\n')
                        print(name, str(length), 'stdOutput: ', stdoutput.strip('\n'), sep=', ', end='***\n')
                        if process.returncode < 0:
                            print('The process calling the slidingWindow program was killed by the system!')
                        print('chenhuilong1\n')
                        failFile.write(name + '\t' + str(length) + '\n')
                        print('***', file=failLogFile)  # 这里同样不用str()也行。都一样。
                        print(process.returncode, name, str(length), sep='\t', file=failLogFile)  # 这里同样不用str()也行。都一样。
                        print('errOutput: ', erroutput.strip('\n'), file=failLogFile)
                        print('stdOutput: ', stdoutput.strip('\n'), end='\n***\n', file=failLogFile)
                        try:
                            os.killpg(process.pid, signal.SIGKILL)
                        except BaseException as e:
                            if str(e) != "[Errno 3] No such process":
                                print(e)

                        removeThreeFiles()
                        # 删除temp.single.input.fasta1.txt，temp.single.input.fasta2.txt， position.txt
                        continue
                    else:  # 子进程状态为0，表明正常执行完毕。根据自己的实验，读写数据量/序列超大的时候，还是可能会出现子进程不被杀死，导致后面读写文件有问题，如position.txt文件出现这行没写完，就停止并写下一行的情况。
                        try:
                            posTable = pd.read_csv('./slidingWindowPosition.txt', encoding='utf-8', sep='\t', header=None)

                            directPositionFile = open(customFolderPath + '/' + name + '_plot_positions_direct.txt', 'w',
                                                      encoding='utf-8')

                            for chl in posTable.values:
                                directPositionFile.write(str(int(chl[0]))+'\t'+str(int(chl[1]))+'\n')
                            directPositionFile.close()
                        except BaseException as e:
                            print(name, str(length), e, e.__traceback__.tb_lineno, sep='***')
                            print('chenhuilong2\n')
                            failFile.write(name + '\t' + str(length) + '\n')
                            print(name, str(length), e, e.__traceback__.tb_lineno, sep='***', file=failLogFile)

                            if os.path.exists(customFolderPath + '/' + name + '_plot_positions_direct.txt'):
                                os.remove(customFolderPath + '/' + name + '_plot_positions_direct.txt')
                            # 如果在，删除directPositionFile，
    
                if '1' in modeList:
                    tempSingleInputFile1 = open('./temp.single.input.fasta1.txt', 'w', encoding='utf-8')
                    tempSingleInputFile1.write(seq)
                    tempSingleInputFile1.close()
                    tempSingleInputFile2 = open('./temp.single.input.fasta2.txt', 'w', encoding='utf-8')
                    tempSingleInputFile2.write(seq)
                    tempSingleInputFile2.close()
                    mode = 1
                    # command = 'cd ./bin/; ./slidingWindow' + ' ' + str(ideSimThr) + ' ' + str(
                    #     windowSize) + ' ' + str(mode) + ' ' + str(aminoAcidMatrix) + ' ' + str(fileType)
                    # os.system(command)
                    command = "./slidingWindow {} {} {} {} {}".format(str(ideSimThr), str(windowSize), str(mode), str(aminoAcidMatrix), str(fileType))
                    process = subprocess.Popen(command, universal_newlines=True, stdin=subprocess.PIPE,
                                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,
                                               preexec_fn=os.setsid)
                    stdoutput, erroutput = process.communicate()  # 这行代码保证调用的子进程结束之后再执行Python脚本中下面的代码。

                    print("mode 1, {}, process.returncode: {}".format(name, process.returncode))
                    if process.returncode:  # 获取进程的返回值。如果进程还没有结束，返回None。根据自我实验，就是等价于Popen.poll()
                        print(name, str(length), 'errOutput: ', erroutput.strip('\n'), sep=', ', end='***\n')
                        print(name, str(length), 'stdOutput: ', stdoutput.strip('\n'), sep=', ', end='***\n')
                        if process.returncode < 0:
                            print('The process calling the slidingWindow program was killed by the system!')
                        print('chenhuilong1\n')
                        failFile.write(name + '\t' + str(length) + '\n')
                        print('***', file=failLogFile)  # 这里同样不用str()也行。都一样。
                        print(process.returncode, name, str(length), sep='\t', file=failLogFile)  # 这里同样不用str()也行。都一样。
                        print('errOutput: ', erroutput.strip('\n'), file=failLogFile)
                        print('stdOutput: ', stdoutput.strip('\n'), end='\n***\n', file=failLogFile)
                        try:
                            os.killpg(process.pid, signal.SIGKILL)
                        except BaseException as e:
                            if str(e) != "[Errno 3] No such process":
                                print(e)

                        removeThreeFiles()
                        # 删除temp.single.input.fasta1.txt，temp.single.input.fasta2.txt， position.txt
                        continue
                    else:
                        try:
                            posTable = pd.read_csv('./slidingWindowPosition.txt', encoding='utf-8', sep='\t', header=None)

                            invertedPositionFile = open(customFolderPath + '/' + name + '_plot_positions_inverted.txt', 'w',
                                                        encoding='utf-8')

                            for chl in posTable.values:
                                invertedPositionFile.write(str(int(chl[0]))+'\t'+str(int(chl[1]))+'\n')
                            invertedPositionFile.close()
                        except BaseException as e:
                            print(name, str(length), e, e.__traceback__.tb_lineno, sep='***')
                            print('chenhuilong2\n')
                            failFile.write(name + '\t' + str(length) + '\n')
                            print(name, str(length), e, e.__traceback__.tb_lineno, sep='***', file=failLogFile)

                            if os.path.exists(customFolderPath + '/' + name + '_plot_positions_inverted.txt'):
                                os.remove(customFolderPath + '/' + name + '_plot_positions_inverted.txt')
                            # 如果在，invertedPositionFile.
            except BaseException as e:
                print(name, str(length), e, e.__traceback__.tb_lineno, sep='***')
                print('chenhuilong3\n')
                failFile.write(name + '\t' + str(length) + '\n')
                print(name, str(length), e, e.__traceback__.tb_lineno, sep='***', file=failLogFile)

            removeThreeFiles()
            # 删除temp.single.input.fasta1.txt，temp.single.input.fasta2.txt， position.txt

        logFile.close()
        failFile.close()
        failLogFile.close()
        seqLengthFileOriginal.close()
    
    for fileName1 in importFileList:
        filePath1 = os.path.join(importFolderPath, fileName1)
        fastaFile = open(filePath1, 'r', encoding='utf-8')
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
        fastaFile.close()
        exportFolderPositionPath = exportFolderPath + '/' + fileName1
        batchExportPosition(chlFasta,exportFolderPositionPath,ideSimThr,windowSize,modeList,aminoAcidMatrix,fileType)

        # 如果结果位置文件的ID在失败的ID中，那删除这个位置文件——做一个检查。
        totalNumber = len(list(chlFasta.keys()))

        folder = '/plot_positions'
        customFailFilePath = exportFolderPositionPath + folder + '/failure_sequences.txt'
        with open(customFailFilePath, 'r', encoding='utf-8') as failFile:
            deduplicationFailNameDict = {}
            for line in failFile:
                print(line.strip('\n'))
                lineList = line.strip('\n').split('\t')
                deduplicationFailNameDict[lineList[0]] = ''
            failNameList = list(deduplicationFailNameDict.keys())
            failNumber = len(failNameList)
        # failFile.close()

        # 'plot_positions' folder
        customFolderPath = exportFolderPositionPath + folder
        customFileList = os.listdir(customFolderPath)

        deduplicationDict = {}
        for fileName in customFileList:
            if '_plot_positions_' in fileName:
                name = re.findall("(.*?)_plot_positions_", fileName)[0]
                deduplicationDict[name] = ''
                if name in failNameList:
                    print("For 'plot_positions' folder, This '{}' is in 'failure_sequences.txt'. it's illogical!".format(name))
                    os.remove(customFolderPath + '/' + fileName)

        successNumber = len(list(deduplicationDict.keys()))

        customFileList = os.listdir(customFolderPath)
        deduplicationAfterDeleteDict = {}
        for fileName in customFileList:
            if '_plot_positions_' in fileName:
                name = re.findall("(.*?)_plot_positions_", fileName)[0]
                deduplicationAfterDeleteDict[name] = ''

        successNumberAfterDelete = len(list(deduplicationAfterDeleteDict.keys()))

        print("For 'plot_positions',")
        print("This fasta file with {} sequences has a total of {} failures and {} successes this time!".format(
            totalNumber, failNumber, successNumber))
        print(
            "If a successful ID appears in 'failure_sequences.txt', the number result after deleting the position file corresponding to this ID is: \n "
            "This fasta file with {} sequences has a total of {} failures and {} successes this time!".format(
                totalNumber, failNumber, successNumberAfterDelete))

        if totalNumber == failNumber + successNumber:
            print(
                "For 'plot_positions', the sum of this number of successes and failures is correct!")
        else:
            state = -1

        # 'positions' folder
        intersectionList = list(set(failNameList) & set(list(deduplicationDict.keys())))
        unionList = list(set(failNameList) | set(list(deduplicationDict.keys())))
        besidesIntersectionList = list(set(unionList) ^ set(list(chlFasta.keys())))
        # if [] == intersectionList and [] == besidesIntersectionList:
        #     print('It is right!')

        if [] == intersectionList and [] == besidesIntersectionList and totalNumber == failNumber + successNumber:
            print(
                "For 'plot_positions', the sum of this number of successes and failures and and the ID is correct!\n")
        else:
            state = -1

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

# created by Huilong Chen, May 24 2022!
# revised by Huilong Chen, May 24 2022!
# revised by Huilong Chen, July 27, 2022! A new CyDotian -> seqAlignToolkit.
# revised by Huilong Chen, August 1, 2022! Optimize.
