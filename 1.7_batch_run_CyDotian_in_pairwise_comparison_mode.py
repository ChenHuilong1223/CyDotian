import os, re, time, argparse, codecs, sys, subprocess, signal
import pandas as pd
import matplotlib.pyplot as plt
sys.stdout = codecs.getwriter("utf-8")(sys.stdout.detach())

parser = argparse.ArgumentParser(description="Batch run CyDotian's modules")
parser.add_argument('-c', '--inputConfigFilePath', metavar='config', required=True, help='Input config file path')
parser.add_argument('-v', '--inputVerticalSequenceFilePath', metavar='vertical', required=True, help='Input vertical sequence file path')
parser.add_argument('-l', '--inputHorizontalSequenceFilePath', metavar='horizontal', required=True, help='Input horizontal sequence file path')
parser.add_argument('-o', '--outputFolderPath', metavar='output', required=True, help='Output folder path')
args = parser.parse_args()

startTime = time.time()
state = 0
try:
    plt.rcParams['font.family'] = ["Times New Roman"]
    fontDict = {"size": 13, "color": "k", 'family': 'Times New Roman'}

    fastaFile1 = open(args.inputVerticalSequenceFilePath, 'r', encoding='utf-8')
    chlVerticalFasta = {}
    gene = seq = ''
    for row in fastaFile1:
        row = row.strip('\n')
        if row.startswith('>'):
            if gene != '' and seq != '':
                chlVerticalFasta[gene] = seq.upper()
            gene = row.replace('>', '')
            seq = ''
        else:
            seq += row
    chlVerticalFasta[gene] = seq.upper()
    fastaFile1.close()

    fastaFile2 = open(args.inputHorizontalSequenceFilePath, 'r', encoding='utf-8')
    chlHorizontalFasta = {}
    gene = seq = ''
    for row in fastaFile2:
        row = row.strip('\n')
        if row.startswith('>'):
            if gene != '' and seq != '':
                chlHorizontalFasta[gene] = seq.upper()
            gene = row.replace('>', '')
            seq = ''
        else:
            seq += row
    chlHorizontalFasta[gene] = seq.upper()
    fastaFile2.close()

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
    fileType,DNA_Matrix,aminoAcidMatrix,modeList,identityThr,similarityThr,repeatLen = 'DNA', 0, 1, ['0', '1'], 0.85, 0.75, 4
    print('=====================')
    for row in configFile:
        if False == row.startswith('#'):
            if '=' in row:
                parameterName = re.findall("(.*?)=", row.strip('\n'))[0]
                matchStr = re.findall("=(.*?)#", row.strip('\n'))[0]
                print(parameterName, matchStr, sep='=')
                parameter = matchStr.replace(' ','')
                if row.startswith('fileType'):
                    if '0' == parameter:
                        fileType = 'DNA'
                    elif '1' == parameter:
                        fileType = 'Amino acid'
                elif row.startswith('DNA_Matrix'):
                    DNA_Matrix = int(parameter)
                elif row.startswith('aminoAcidMatrix'):
                    aminoAcidMatrix = int(parameter)
                elif row.startswith('mode'):
                    modeList = parameter.split(',')
                elif row.startswith('identityThr'):
                    identityThr = float(parameter)
                elif row.startswith('similarityThr'):
                    similarityThr = float(parameter)
                elif row.startswith('repeatLen'):
                    repeatLen = int(parameter)
    configFile.close()
    print('=====================\n')

    def batchExportPosition(chlVerticalFasta,chlHorizontalFasta, exportFolderPath, fileType, DNA_Matrix, aminoAcidMatrix, modeList,
                            identityThr, similarityThr, repeatLen):
        folderOriginal = '/positions_original'
        customFolderPathOriginal = exportFolderPath + folderOriginal
        isExists = os.path.exists(customFolderPathOriginal)
        if isExists:
            pass
        else:
            os.makedirs(customFolderPathOriginal, mode=0o777)

        logFile = open(customFolderPathOriginal + '/log.txt', 'w', encoding='utf-8')
        if 'DNA' == fileType:
            DNA_MatrixList = ['BLAST matrix', 'Transition-transversion matrix']
            logFile.write('# DNA' + '\t' + DNA_MatrixList[DNA_Matrix] + '\n')
            logFile.write(
                '# parameter: Identity >= ' + str(identityThr * 100) + '% ' + 'RepeatLen >= ' + str(repeatLen) + '\n')
            logFile.write(
                '# header of positions.txt: Start1  End1    Start2  End2    Length   Identity    Mismatch    Score' + '\n')
        elif 'Amino acid' == fileType:
            aminoAcidMatrixList = ['BLOSUM45', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90', 'PAM30', 'PAM70', 'PAM250']
            logFile.write('# Amino acid' + '\t' + aminoAcidMatrixList[aminoAcidMatrix] + '\n')
            logFile.write(
                '# parameter: Similarity >= ' + str(similarityThr * 100) + '% ' + 'RepeatLen >= ' + str(
                    repeatLen) + '\n')
            logFile.write(
                '# header of positions.txt: Start1  End1    Start2  End2    Length   Identity    Mismatch    '
                'Similarity  Score' + '\n')
        logFile.write('# Sequence names containing repeats larger than the RepeatLen threshold are recorded in file '
                      'sequence_names_direct.txt or sequence_names_inverted.txt.'
                      ' The names and lengths of sequences that failed to be analysed are recorded in file failure_sequences.txt.\n')
        logFile.write('Time: ' + time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))

        failFile = open(customFolderPathOriginal + '/failure_sequences.txt', 'w', encoding='utf-8')
        failLogFile = open(customFolderPathOriginal + '/failure_sequences_error_log.txt', 'w', encoding='utf-8')

        if '0' in modeList:
            directPositionBigLengthLogFile = open(customFolderPathOriginal + '/sequence_names_direct.txt', 'w', encoding='utf-8')
        if '1' in modeList:
            invertedPositionBigLengthLogFile = open(customFolderPathOriginal + '/sequence_names_inverted.txt', 'w', encoding='utf-8')

        seqLength1FileOriginal = open(customFolderPathOriginal + '/all_vertical_sequences_length.txt', 'w', encoding='utf-8')
        seqLength2FileOriginal = open(customFolderPathOriginal + '/all_horizontal_sequences_length.txt', 'w',
                                      encoding='utf-8')
        for name1, seq1 in chlVerticalFasta.items():
            seqLength1FileOriginal.write(name1 + '\t' + str(len(seq1)) + '\n')

        for name2, seq2 in chlHorizontalFasta.items():
            seqLength2FileOriginal.write(name2 + '\t' + str(len(seq2)) + '\n')

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

        for name1, seq1 in chlVerticalFasta.items():
            length1 = len(seq1)
            for name2, seq2 in chlHorizontalFasta.items():
                length2 = len(seq2)
                try:
                    name1vs2 = name1 + '_VS_' + name2
                    if '0' in modeList:
                        if 'DNA' == fileType:
                            tempSingleInputFile1 = open('./temp.single.input.fasta1.txt', 'w', encoding='utf-8')
                            tempSingleInputFile1.write(seq1)
                            tempSingleInputFile1.close()
                            tempSingleInputFile2 = open('./temp.single.input.fasta2.txt', 'w', encoding='utf-8')
                            tempSingleInputFile2.write(seq2)
                            tempSingleInputFile2.close()
                            mode = 0
                            exactMatch = -1

                            command = "./bpRepeatScan {} {} {} {}".format(str(identityThr), str(mode), str(DNA_Matrix), str(exactMatch))
                            process = subprocess.Popen(command, universal_newlines=True, stdin=subprocess.PIPE,
                                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,
                                                       preexec_fn=os.setsid)
                            stdoutput, erroutput = process.communicate()

                            print("mode 0, {}, process.returncode: {}".format(name1vs2, process.returncode))
                            if process.returncode:
                                print(name1vs2, length1, length2, 'errOutput: ', erroutput.strip('\n'), sep=', ', end='***\n')
                                print(name1vs2, length1, length2, 'stdOutput: ', stdoutput.strip('\n'), sep=', ', end='***\n')
                                if process.returncode < 0:
                                    print('The process calling the bpRepeatScan program was killed by the system!')
                                print('chenhuilong1\n')
                                failFile.write(name1vs2 + '\t' + str(length1) + '\t' + str(length2) + '\n')
                                print('***', file=failLogFile)
                                print(process.returncode, name1vs2, length1, length2, sep='\t', end='***\n', file=failLogFile)
                                print('errOutput: ', erroutput.strip('\n'), file=failLogFile)
                                print('stdOutput: ', stdoutput.strip('\n'), end='\n***\n', file=failLogFile)
                                try:
                                    os.killpg(process.pid, signal.SIGKILL)
                                except BaseException as e:
                                    if str(e) != "[Errno 3] No such process":
                                        print(e)

                                removeThreeFiles()
                                continue
                            else:
                                try:
                                    posTable = pd.read_csv('./position.txt', encoding='utf-8', sep='\t', header=None)
                                    posTable = posTable.sort_values(by=2, ascending=False).reset_index(drop=True)

                                    if posTable[2][1] >= repeatLen:
                                        directPositionBigLengthLogFile.write(name1vs2 + '\n')

                                    directPositionFileOriginal = open(
                                        customFolderPathOriginal + '/' + name1vs2 + '_positions_direct.txt', 'w',
                                        encoding='utf-8')

                                    for chl in posTable.values:
                                        if chl[2] >= repeatLen:
                                            directPositionFileOriginal.write(
                                                str(int(chl[0])) + '\t' + str(int(chl[0] + chl[2] - 1)) + '\t' + str(int(chl[1])) + '\t'
                                                + str(int(chl[1] + chl[2] - 1)) + '\t' + str(int(chl[2])) + '\t' + str(
                                                    '{:6f}'.format(chl[3]))
                                                + '\t' + str(int(chl[4])) + '\t' + str(int(chl[5])) + '\n')
                                    directPositionFileOriginal.close()
                                except BaseException as e:
                                    print(name1vs2, length1, length2, e, e.__traceback__.tb_lineno, sep='***')
                                    print('chenhuilong2\n')
                                    failFile.write(name1vs2 + '\t' + str(length1) + '\t' + str(length2) + '\n')
                                    print(name1vs2, length1, length2, e, e.__traceback__.tb_lineno, sep='***', file=failLogFile)

                                    if os.path.exists(customFolderPathOriginal + '/' + name1vs2 + '_positions_direct.txt'):
                                        os.remove(customFolderPathOriginal + '/' + name1vs2 + '_positions_direct.txt')
                        elif 'Amino acid' == fileType:
                            tempSingleInputFile1 = open('./temp.single.input.fasta1.txt', 'w', encoding='utf-8')
                            tempSingleInputFile1.write(seq1)
                            tempSingleInputFile1.close()
                            tempSingleInputFile2 = open('./temp.single.input.fasta2.txt', 'w', encoding='utf-8')
                            tempSingleInputFile2.write(seq2)
                            tempSingleInputFile2.close()
                            mode = 0
                            exactMatch = -1

                            command = "./aaRepeatScan {} {} {} {}".format(str(similarityThr), str(mode), str(aminoAcidMatrix), str(exactMatch))
                            process = subprocess.Popen(command, universal_newlines=True, stdin=subprocess.PIPE,
                                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,
                                                       preexec_fn=os.setsid)
                            stdoutput, erroutput = process.communicate()

                            print("mode 0, {}, process.returncode: {}".format(name1vs2, process.returncode))
                            if process.returncode:
                                print(name1vs2, length1, length2, 'errOutput: ', erroutput.strip('\n'), sep=', ', end='***\n')
                                print(name1vs2, length1, length2, 'stdOutput: ', stdoutput.strip('\n'), sep=', ', end='***\n')
                                if process.returncode < 0:
                                    print('The process calling the aaRepeatScan program was killed by the system!')
                                print('chenhuilong1\n')
                                failFile.write(name1vs2 + '\t' + str(length1) + '\t' + str(length2) + '\n')
                                print('***', file=failLogFile)
                                print(process.returncode, name1vs2, length1, length2, sep='\t', end='***\n',
                                      file=failLogFile)
                                print('errOutput: ', erroutput.strip('\n'), file=failLogFile)
                                print('stdOutput: ', stdoutput.strip('\n'), end='\n***\n', file=failLogFile)
                                try:
                                    os.killpg(process.pid, signal.SIGKILL)
                                except BaseException as e:
                                    if str(e) != "[Errno 3] No such process":
                                        print(e)

                                removeThreeFiles()
                                continue
                            else:
                                try:
                                    posTable = pd.read_csv('./position.txt', encoding='utf-8', sep='\t', header=None)
                                    posTable = posTable.sort_values(by=2, ascending=False).reset_index(drop=True)

                                    if posTable[2][1] >= repeatLen:
                                        directPositionBigLengthLogFile.write(name1vs2 + '\n')

                                    directPositionFileOriginal = open(
                                        customFolderPathOriginal + '/' + name1vs2 + '_positions_direct.txt', 'w',
                                        encoding='utf-8')

                                    for chl in posTable.values:
                                        if chl[2] >= repeatLen:
                                            directPositionFileOriginal.write(
                                                str(int(chl[0])) + '\t' + str(int(chl[0] + chl[2] - 1)) + '\t' + str(int(chl[1])) + '\t'
                                                + str(int(chl[1] + chl[2] - 1)) + '\t' + str(int(chl[2])) + '\t' + str(
                                                    '{:6f}'.format(chl[3]))
                                                + '\t' + str(int(chl[4])) + '\t' + str('{:6f}'.format(chl[5])) + '\t' + str(
                                                    int(chl[6])) + '\n')
                                    directPositionFileOriginal.close()
                                except BaseException as e:
                                    print(name1vs2, length1, length2, e, e.__traceback__.tb_lineno, sep='***')
                                    print('chenhuilong2\n')
                                    failFile.write(name1vs2 + '\t' + str(length1) + '\t' + str(length2) + '\n')
                                    print(name1vs2, length1, length2, e, e.__traceback__.tb_lineno, sep='***', file=failLogFile)

                                    if os.path.exists(customFolderPathOriginal + '/' + name1vs2 + '_positions_direct.txt'):
                                        os.remove(customFolderPathOriginal + '/' + name1vs2 + '_positions_direct.txt')

                    if '1' in modeList:
                        if 'DNA' == fileType:
                            tempSingleInputFile1 = open('./temp.single.input.fasta1.txt', 'w', encoding='utf-8')
                            tempSingleInputFile1.write(seq1)
                            tempSingleInputFile1.close()
                            tempSingleInputFile2 = open('./temp.single.input.fasta2.txt', 'w', encoding='utf-8')
                            tempSingleInputFile2.write(seq2)
                            tempSingleInputFile2.close()
                            mode = 1
                            exactMatch = -1

                            command = "./bpRepeatScan {} {} {} {}".format(str(identityThr), str(mode), str(DNA_Matrix), str(exactMatch))
                            process = subprocess.Popen(command, universal_newlines=True, stdin=subprocess.PIPE,
                                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,
                                                       preexec_fn=os.setsid)
                            stdoutput, erroutput = process.communicate()

                            print("mode 1, {}, process.returncode: {}".format(name1vs2, process.returncode))
                            if process.returncode:
                                print(name1vs2, length1, length2, 'errOutput: ', erroutput.strip('\n'), sep=', ', end='***\n')
                                print(name1vs2, length1, length2, 'stdOutput: ', stdoutput.strip('\n'), sep=', ', end='***\n')
                                if process.returncode < 0:
                                    print('The process calling the bpRepeatScan program was killed by the system!')
                                print('chenhuilong1\n')
                                failFile.write(name1vs2 + '\t' + str(length1) + '\t' + str(length2) + '\n')
                                print('***', file=failLogFile)
                                print(process.returncode, name1vs2, length1, length2, sep='\t', end='***\n',
                                      file=failLogFile)
                                print('errOutput: ', erroutput.strip('\n'), file=failLogFile)
                                print('stdOutput: ', stdoutput.strip('\n'), end='\n***\n', file=failLogFile)
                                try:
                                    os.killpg(process.pid, signal.SIGKILL)
                                except BaseException as e:
                                    if str(e) != "[Errno 3] No such process":
                                        print(e)

                                removeThreeFiles()
                                continue
                            else:
                                try:
                                    posTable = pd.read_csv('./position.txt', encoding='utf-8', sep='\t', header=None)
                                    posTable = posTable.sort_values(by=2, ascending=False).reset_index(drop=True)

                                    if posTable[2][0] >= repeatLen:
                                        invertedPositionBigLengthLogFile.write(name1vs2 + '\n')

                                    invertedPositionFileOriginal = open(
                                        customFolderPathOriginal + '/' + name1vs2 + '_positions_inverted.txt',
                                        'w', encoding='utf-8')

                                    for chl in posTable.values:
                                        if chl[2] >= repeatLen:
                                            invertedPositionFileOriginal.write(
                                                str(int(chl[0])) + '\t' + str(int(chl[0] + chl[2] - 1)) + '\t' + str(int(chl[1])) + '\t'
                                                + str(int(chl[1] - chl[2] + 1)) + '\t' + str(int(chl[2])) + '\t' + str(
                                                    '{:6f}'.format(chl[3]))
                                                + '\t' + str(int(chl[4])) + '\t' + str(int(chl[5])) + '\n')
                                    invertedPositionFileOriginal.close()
                                except BaseException as e:
                                    print(name1vs2, length1, length2, e, e.__traceback__.tb_lineno, sep='***')
                                    print('chenhuilong2\n')
                                    failFile.write(name1vs2 + '\t' + str(length1) + '\t' + str(length2) + '\n')
                                    print(name1vs2, length1, length2, e, e.__traceback__.tb_lineno, sep='***', file=failLogFile)

                                    if os.path.exists(customFolderPathOriginal + '/' + name1vs2 + '_positions_inverted.txt'):
                                        os.remove(customFolderPathOriginal + '/' + name1vs2 + '_positions_inverted.txt')
                        elif 'Amino acid' == fileType:
                            tempSingleInputFile1 = open('./temp.single.input.fasta1.txt', 'w', encoding='utf-8')
                            tempSingleInputFile1.write(seq1)
                            tempSingleInputFile1.close()
                            tempSingleInputFile2 = open('./temp.single.input.fasta2.txt', 'w', encoding='utf-8')
                            tempSingleInputFile2.write(seq2)
                            tempSingleInputFile2.close()
                            mode = 1
                            exactMatch = -1

                            command = "./aaRepeatScan {} {} {} {}".format(str(similarityThr), str(mode), str(aminoAcidMatrix), str(exactMatch))
                            process = subprocess.Popen(command, universal_newlines=True, stdin=subprocess.PIPE,
                                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,
                                                       preexec_fn=os.setsid)
                            stdoutput, erroutput = process.communicate()

                            print("mode 1, {}, process.returncode: {}".format(name1vs2, process.returncode))
                            if process.returncode:
                                print(name1vs2, length1, length2, 'errOutput: ', erroutput.strip('\n'), sep=', ', end='***\n')
                                print(name1vs2, length1, length2, 'stdOutput: ', stdoutput.strip('\n'), sep=', ', end='***\n')
                                if process.returncode < 0:
                                    print('The process calling the aaRepeatScan program was killed by the system!')
                                print('chenhuilong1\n')
                                failFile.write(name1vs2 + '\t' + str(length1) + '\t' + str(length2) + '\n')
                                print('***', file=failLogFile)
                                print(process.returncode, name1vs2, length1, length2, sep='\t', end='***\n',
                                      file=failLogFile)
                                print('errOutput: ', erroutput.strip('\n'), file=failLogFile)
                                print('stdOutput: ', stdoutput.strip('\n'), end='\n***\n', file=failLogFile)
                                try:
                                    os.killpg(process.pid, signal.SIGKILL)
                                except BaseException as e:
                                    if str(e) != "[Errno 3] No such process":
                                        print(e)

                                removeThreeFiles()
                                continue
                            else:
                                try:
                                    posTable = pd.read_csv('./position.txt', encoding='utf-8', sep='\t', header=None)
                                    posTable = posTable.sort_values(by=2, ascending=False).reset_index(drop=True)

                                    if posTable[2][0] >= repeatLen:
                                        invertedPositionBigLengthLogFile.write(name1vs2 + '\n')

                                    invertedPositionFileOriginal = open(
                                        customFolderPathOriginal + '/' + name1vs2 + '_positions_inverted.txt',
                                        'w', encoding='utf-8')

                                    for chl in posTable.values:
                                        if chl[2] >= repeatLen:
                                            invertedPositionFileOriginal.write(
                                                str(int(chl[0])) + '\t' + str(int(chl[0] + chl[2] - 1)) + '\t' + str(int(chl[1])) + '\t'
                                                + str(int(chl[1] - chl[2] + 1)) + '\t' + str(int(chl[2])) + '\t' + str(
                                                    '{:6f}'.format(chl[3]))
                                                + '\t' + str(int(chl[4])) + '\t' + str('{:6f}'.format(chl[5])) + '\t' + str(
                                                    int(chl[6])) + '\n')
                                    invertedPositionFileOriginal.close()
                                except BaseException as e:
                                    print(name1vs2, length1, length2, e, e.__traceback__.tb_lineno, sep='***')
                                    print('chenhuilong2\n')
                                    failFile.write(name1vs2 + '\t' + str(length1) + '\t' + str(length2) + '\n')
                                    print(name1vs2, length1, length2, e, e.__traceback__.tb_lineno, sep='***', file=failLogFile)

                                    if os.path.exists(customFolderPathOriginal + '/' + name1vs2 + '_positions_inverted.txt'):
                                        os.remove(customFolderPathOriginal + '/' + name1vs2 + '_positions_inverted.txt')
                except BaseException as e:
                    print(name1vs2, length1, length2, e, e.__traceback__.tb_lineno, sep='***')
                    print('chenhuilong3\n')
                    failFile.write(name1vs2 + '\t' + str(length1) + '\t' + str(length2) + '\n')
                    print(name1vs2, length1, length2, e, e.__traceback__.tb_lineno, sep='***', file=failLogFile)

                removeThreeFiles()

        logFile.close()
        failFile.close()
        failLogFile.close()
        if '0' in modeList:
            directPositionBigLengthLogFile.close()
        if '1' in modeList:
            invertedPositionBigLengthLogFile.close()
        seqLength1FileOriginal.close()
        seqLength2FileOriginal.close()

    batchExportPosition(chlVerticalFasta, chlHorizontalFasta, exportFolderPath, fileType, DNA_Matrix, aminoAcidMatrix, modeList,
                        identityThr, similarityThr, repeatLen)

    def directDotplot(filePath):
        directPosFile = open(filePath, 'r', encoding='utf-8')
        x, y = [], []
        for line in directPosFile:
            lineList = line.strip('\n').split('\t')
            for i in range(int(lineList[0]), int(lineList[0]) + int(lineList[4])):
                y.append(i)
            for j in range(int(lineList[2]), int(lineList[2]) + int(lineList[4])):
                x.append(j)
        directPosFile.close()
        return x, y
    def invertedDotplot(filePath):
        invertedPosFile = open(filePath, 'r', encoding='utf-8')
        x, y = [], []
        for line in invertedPosFile:
            lineList = line.strip('\n').split('\t')
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
            lineList = line.strip('\n').rsplit('\t', 1)
            seqLengthDic[lineList[0]] = lineList[1]
        seqLengthFile.close()
        return seqLengthDic

    positionsOriginalFolder = exportFolderPath + '/positions_original'
    positionsOriginalFileList = os.listdir(positionsOriginalFolder)
    positionsOriginalFolderPath = os.path.abspath(positionsOriginalFolder)

    exportFolderDotplot = exportFolderPath + '/positions_original_dotplots'
    mkdir(exportFolderDotplot)
    exportFolderDotplotPath = os.path.abspath(exportFolderDotplot)

    logFile = open(exportFolderDotplotPath + '/log.txt', 'w', encoding='utf-8')
    logFile.write('# The names of files that failed to be analysed are recorded in file failure_files.txt.\n')
    logFile.write('Time: ' + time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))

    failFile = open(exportFolderDotplotPath + '/failure_files.txt', 'w', encoding='utf-8')

    seqLength1FilePath = os.path.join(positionsOriginalFolderPath, 'all_vertical_sequences_length.txt')
    seqLength1Dic = creatSeqLengthDic(seqLength1FilePath)
    seqLength2FilePath = os.path.join(positionsOriginalFolderPath, 'all_horizontal_sequences_length.txt')
    seqLength2Dic = creatSeqLengthDic(seqLength2FilePath)

    for fileName in positionsOriginalFileList:
        try:
            filePath = os.path.join(positionsOriginalFolderPath, fileName)
            if fileName.endswith('_positions_direct.txt'):
                name = re.findall("(.*?)_positions_direct.txt", fileName)[0]
                name1, name2 = name.split('_VS_')[0], name.split('_VS_')[1]
                length1 = int(seqLength1Dic[name1])
                length2 = int(seqLength2Dic[name2])
                fig, ax = plt.subplots()
                ax.xaxis.set_ticks_position('top')
                ax.invert_yaxis()
                ax.set_aspect(1)
                x, y = directDotplot(filePath)
                ax.scatter(x, y, marker='o', c='k', s=0.5, linewidths=0, edgecolors=None, alpha=1)
                ax.set_xlim(0, length2 + 1)
                ax.set_ylim(length1 + 1, 0)
                ax.set_xlabel(name2, fontdict=fontDict)
                ax.xaxis.set_label_position('top')
                ax.set_ylabel(name1, fontdict=fontDict)
                fig.savefig(exportFolderDotplotPath + '/' + name + '_direct.png', format='png', dpi=500)
                fig.savefig(exportFolderDotplotPath + '/' + name + '_direct.pdf', format='pdf')
                plt.clf()
            elif fileName.endswith('_positions_inverted.txt'):
                name = re.findall("(.*?)_positions_inverted.txt", fileName)[0]
                name1, name2 = name.split('_VS_')[0], name.split('_VS_')[1]
                length1 = int(seqLength1Dic[name1])
                length2 = int(seqLength2Dic[name2])
                fig, ax = plt.subplots()
                ax.xaxis.set_ticks_position('top')
                ax.invert_yaxis()
                ax.set_aspect(1)
                x, y = invertedDotplot(filePath)
                ax.scatter(x, y, marker='o', c='k', s=0.5, linewidths=0, edgecolors=None, alpha=1)
                ax.set_xlim(0, length2 + 1)
                ax.set_ylim(length1 + 1, 0)
                ax.set_xlabel(name2, fontdict=fontDict)
                ax.xaxis.set_label_position('top')
                ax.set_ylabel(name1, fontdict=fontDict)
                fig.savefig(exportFolderDotplotPath + '/' + name + '_inverted.png', format='png', dpi=500)
                fig.savefig(exportFolderDotplotPath + '/' + name + '_inverted.pdf', format='pdf')
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

# created by Huilong Chen, May 27, 2022!
# revised by Huilong Chen, May 27, 2022!
# revised by Huilong Chen, July 26, 2022! A new CyDotian.
# revised by Huilong Chen, August 1, 2022! Optimize.
