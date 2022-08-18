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

    def batchExportPosition(chlFasta,exportFolderPositionPath,fileType,DNA_Matrix,aminoAcidMatrix,modeList,identityThr,similarityThr,repeatLen):
        folder = '/positions'
        customFolderPath = exportFolderPositionPath + folder
        isExists = os.path.exists(customFolderPath)
        if isExists:
            pass
        else:
            os.makedirs(customFolderPath, mode=0o777)

        folderOriginal = '/positions_original'
        customFolderPathOriginal = exportFolderPositionPath + folderOriginal
        isExists = os.path.exists(customFolderPathOriginal)
        if isExists:
            pass
        else:
            os.makedirs(customFolderPathOriginal, mode=0o777)

        logFile = open(customFolderPath + '/log.txt', 'w',encoding='utf-8')
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
                      'sequence_names_direct.txt, sequence_names_inverted.txt or sequence_names_reverse_complement.txt.'
                      ' The names and lengths of sequences that failed to be analysed are recorded in file failure_sequences.txt.\n')
        logFile.write('Time: ' + time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    
        failFile = open(customFolderPath + '/failure_sequences.txt', 'w',encoding='utf-8')
        failLogFile = open(customFolderPath + '/failure_sequences_error_log.txt', 'w',encoding='utf-8')

        if '0' in modeList:
            directPositionBigLengthLogFile = open(customFolderPath + '/sequence_names_direct.txt', 'w',encoding='utf-8')
        if '1' in modeList:
            invertedPositionBigLengthLogFile = open(customFolderPath + '/sequence_names_inverted.txt','w',encoding='utf-8')
        if '2' in modeList:
            reverseComplementPositionBigLengthLogFile = open(
                customFolderPath + '/sequence_names_reverse_complement.txt', 'w', encoding='utf-8')

        seqLengthFileOriginal = open(customFolderPathOriginal + '/all_sequences_length.txt', 'w', encoding='utf-8')

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
                    if 'DNA' == fileType:
                        tempSingleInputFile1 = open('./temp.single.input.fasta1.txt','w',encoding='utf-8')
                        tempSingleInputFile1.write(seq)
                        tempSingleInputFile1.close()
                        tempSingleInputFile2 = open('./temp.single.input.fasta2.txt','w',encoding='utf-8')
                        tempSingleInputFile2.write(seq)
                        tempSingleInputFile2.close()
                        mode = 0
                        exactMatch = -1

                        command = "./bpRepeatScan {} {} {} {}".format(str(identityThr), str(mode), str(DNA_Matrix), str(exactMatch))
                        process = subprocess.Popen(command, universal_newlines=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, preexec_fn=os.setsid)
                        stdoutput, erroutput = process.communicate()

                        print("mode 0, {}, process.returncode: {}".format(name, process.returncode))
                        if process.returncode:
                            print(name, str(length), 'errOutput: ', erroutput.strip('\n'), sep=', ', end='***\n')
                            print(name, str(length), 'stdOutput: ', stdoutput.strip('\n'), sep=', ', end='***\n')
                            if process.returncode < 0:
                                print('The process calling the bpRepeatScan program was killed by the system!')
                            print('chenhuilong1\n')
                            failFile.write(name + '\t' + str(length) + '\n')
                            print('***', file=failLogFile)
                            print(process.returncode, name, str(length), sep='\t', file=failLogFile)
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
                                    directPositionBigLengthLogFile.write(name + '\n')

                                directPositionFile = open(customFolderPath + '/' + name + '_positions_direct.txt', 'w',
                                                          encoding='utf-8')
                                directPositionFileOriginal = open(
                                    customFolderPathOriginal + '/' + name + '_positions_direct.txt', 'w',
                                    encoding='utf-8')

                                for chl in posTable.values:
                                    if chl[1] > chl[0] and chl[2] >= repeatLen:
                                        directPositionFile.write(
                                            str(int(chl[0])) + '\t' + str(int(chl[0] + chl[2] - 1)) + '\t' + str(int(chl[1])) + '\t'
                                            + str(int(chl[1] + chl[2] - 1)) + '\t' + str(int(chl[2])) + '\t' + str('{:6f}'.format(chl[3]))
                                            + '\t' + str(int(chl[4])) + '\t' + str(int(chl[5])) + '\n')
                                    if chl[2] >= repeatLen:
                                        directPositionFileOriginal.write(
                                            str(int(chl[0])) + '\t' + str(int(chl[0] + chl[2] - 1)) + '\t' + str(int(chl[1])) + '\t'
                                            + str(int(chl[1] + chl[2] - 1)) + '\t' + str(int(chl[2])) + '\t' + str('{:6f}'.format(chl[3]))
                                            + '\t' + str(int(chl[4])) + '\t' + str(int(chl[5])) + '\n')
                                directPositionFile.close()
                                directPositionFileOriginal.close()
                            except BaseException as e:
                                print(name, str(length), e, e.__traceback__.tb_lineno, sep='***')
                                print('chenhuilong2\n')
                                failFile.write(name + '\t' + str(length) + '\n')
                                print(name, str(length), e, e.__traceback__.tb_lineno, sep='***', file=failLogFile)

                                if os.path.exists(customFolderPath + '/' + name + '_positions_direct.txt'):
                                    os.remove(customFolderPath + '/' + name + '_positions_direct.txt')
                                if os.path.exists(customFolderPathOriginal + '/' + name + '_positions_direct.txt'):
                                    os.remove(customFolderPathOriginal + '/' + name + '_positions_direct.txt')
                    elif 'Amino acid' == fileType:
                        tempSingleInputFile1 = open('./temp.single.input.fasta1.txt', 'w', encoding='utf-8')
                        tempSingleInputFile1.write(seq)
                        tempSingleInputFile1.close()
                        tempSingleInputFile2 = open('./temp.single.input.fasta2.txt', 'w', encoding='utf-8')
                        tempSingleInputFile2.write(seq)
                        tempSingleInputFile2.close()
                        mode = 0
                        exactMatch = -1

                        command = "./aaRepeatScan {} {} {} {}".format(str(similarityThr), str(mode), str(aminoAcidMatrix), str(exactMatch))
                        process = subprocess.Popen(command, universal_newlines=True, stdin=subprocess.PIPE,
                                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,
                                                   preexec_fn=os.setsid)
                        stdoutput, erroutput = process.communicate()

                        print("mode 0, {}, process.returncode: {}".format(name, process.returncode))
                        if process.returncode:
                            print(name, str(length), 'errOutput: ', erroutput.strip('\n'), sep=', ', end='***\n')
                            print(name, str(length), 'stdOutput: ', stdoutput.strip('\n'), sep=', ', end='***\n')
                            if process.returncode < 0:
                                print('The process calling the aaRepeatScan program was killed by the system!')
                            print('chenhuilong1\n')
                            failFile.write(name + '\t' + str(length) + '\n')
                            print('***', file=failLogFile)
                            print(process.returncode, name, str(length), sep='\t', file=failLogFile)
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
                                    directPositionBigLengthLogFile.write(name + '\n')

                                directPositionFile = open(customFolderPath + '/' + name + '_positions_direct.txt', 'w',
                                                          encoding='utf-8')
                                directPositionFileOriginal = open(
                                    customFolderPathOriginal + '/' + name + '_positions_direct.txt', 'w',
                                    encoding='utf-8')

                                for chl in posTable.values:
                                    if chl[1] > chl[0] and chl[2] >= repeatLen:
                                        directPositionFile.write(
                                            str(int(chl[0])) + '\t' + str(int(chl[0] + chl[2] - 1)) + '\t' + str(int(chl[1])) + '\t'
                                            + str(int(chl[1] + chl[2] - 1)) + '\t' + str(int(chl[2])) + '\t' + str('{:6f}'.format(chl[3]))
                                            + '\t' + str(int(chl[4])) + '\t' + str('{:6f}'.format(chl[5])) + '\t' + str(int(chl[6])) + '\n')
                                    if chl[2] >= repeatLen:
                                        directPositionFileOriginal.write(
                                            str(int(chl[0])) + '\t' + str(int(chl[0] + chl[2] - 1)) + '\t' + str(int(chl[1])) + '\t'
                                            + str(int(chl[1] + chl[2] - 1)) + '\t' + str(int(chl[2])) + '\t' + str(
                                                '{:6f}'.format(chl[3]))
                                            + '\t' + str(int(chl[4])) + '\t' + str('{:6f}'.format(chl[5])) + '\t' + str(
                                                int(chl[6])) + '\n')
                                directPositionFile.close()
                                directPositionFileOriginal.close()
                            except BaseException as e:
                                print(name, str(length), e, e.__traceback__.tb_lineno, sep='***')
                                print('chenhuilong2\n')
                                failFile.write(name + '\t' + str(length) + '\n')
                                print(name, str(length), e, e.__traceback__.tb_lineno, sep='***', file=failLogFile)

                                if os.path.exists(customFolderPath + '/' + name + '_positions_direct.txt'):
                                    os.remove(customFolderPath + '/' + name + '_positions_direct.txt')
                                if os.path.exists(customFolderPathOriginal + '/' + name + '_positions_direct.txt'):
                                    os.remove(customFolderPathOriginal + '/' + name + '_positions_direct.txt')

                if '1' in modeList:
                    if 'DNA' == fileType:
                        tempSingleInputFile1 = open('./temp.single.input.fasta1.txt', 'w', encoding='utf-8')
                        tempSingleInputFile1.write(seq)
                        tempSingleInputFile1.close()
                        tempSingleInputFile2 = open('./temp.single.input.fasta2.txt', 'w', encoding='utf-8')
                        tempSingleInputFile2.write(seq)
                        tempSingleInputFile2.close()
                        mode = 1
                        exactMatch = -1

                        command = "./bpRepeatScan {} {} {} {}".format(str(identityThr), str(mode), str(DNA_Matrix), str(exactMatch))
                        process = subprocess.Popen(command, universal_newlines=True, stdin=subprocess.PIPE,
                                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,
                                                   preexec_fn=os.setsid)
                        stdoutput, erroutput = process.communicate()

                        print("mode 1, {}, process.returncode: {}".format(name, process.returncode))
                        if process.returncode:
                            print(name, str(length), 'errOutput: ', erroutput.strip('\n'), sep=', ', end='***\n')
                            print(name, str(length), 'stdOutput: ', stdoutput.strip('\n'), sep=', ', end='***\n')
                            if process.returncode < 0:
                                print('The process calling the bpRepeatScan program was killed by the system!')
                            print('chenhuilong1\n')
                            failFile.write(name + '\t' + str(length) + '\n')
                            print('***', file=failLogFile)
                            print(process.returncode, name, str(length), sep='\t', file=failLogFile)
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
                                    invertedPositionBigLengthLogFile.write(name + '\n')

                                invertedPositionFile = open(customFolderPath + '/' + name + '_positions_inverted.txt', 'w',
                                                            encoding='utf-8')
                                invertedPositionFileOriginal = open(
                                    customFolderPathOriginal + '/' + name + '_positions_inverted.txt', 'w', encoding='utf-8')

                                for chl in posTable.values:
                                    if chl[1] > chl[0] and chl[2] >= repeatLen:
                                        invertedPositionFile.write(
                                            str(int(chl[0])) + '\t' + str(int(chl[0] + chl[2] - 1)) + '\t' + str(int(chl[1])) + '\t'
                                            + str(int(chl[1] - chl[2] + 1)) + '\t' + str(int(chl[2])) + '\t' + str(
                                                '{:6f}'.format(chl[3])) + '\t' + str(int(chl[4])) + '\t' + str(int(chl[5])) + '\n')
                                    if chl[2] >= repeatLen:
                                        invertedPositionFileOriginal.write(
                                            str(int(chl[0])) + '\t' + str(int(chl[0] + chl[2] - 1)) + '\t' + str(int(chl[1])) + '\t'
                                            + str(int(chl[1] - chl[2] + 1)) + '\t' + str(int(chl[2])) + '\t' + str(
                                                '{:6f}'.format(chl[3])) + '\t' + str(int(chl[4])) + '\t' + str(int(chl[5])) + '\n')
                                invertedPositionFile.close()
                                invertedPositionFileOriginal.close()
                            except BaseException as e:
                                print(name, str(length), e, e.__traceback__.tb_lineno, sep='***')
                                print('chenhuilong2\n')
                                failFile.write(name + '\t' + str(length) + '\n')
                                print(name, str(length), e, e.__traceback__.tb_lineno, sep='***', file=failLogFile)

                                if os.path.exists(customFolderPath + '/' + name + '_positions_inverted.txt'):
                                    os.remove(customFolderPath + '/' + name + '_positions_inverted.txt')
                                if os.path.exists(customFolderPathOriginal + '/' + name + '_positions_inverted.txt'):
                                    os.remove(customFolderPathOriginal + '/' + name + '_positions_inverted.txt')
                    elif 'Amino acid' == fileType:
                        tempSingleInputFile1 = open('./temp.single.input.fasta1.txt', 'w', encoding='utf-8')
                        tempSingleInputFile1.write(seq)
                        tempSingleInputFile1.close()
                        tempSingleInputFile2 = open('./temp.single.input.fasta2.txt', 'w', encoding='utf-8')
                        tempSingleInputFile2.write(seq)
                        tempSingleInputFile2.close()
                        mode = 1
                        exactMatch = -1

                        command = "./aaRepeatScan {} {} {} {}".format(str(similarityThr), str(mode), str(aminoAcidMatrix), str(exactMatch))
                        process = subprocess.Popen(command, universal_newlines=True, stdin=subprocess.PIPE,
                                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,
                                                   preexec_fn=os.setsid)
                        stdoutput, erroutput = process.communicate()

                        print("mode 1, {}, process.returncode: {}".format(name, process.returncode))
                        if process.returncode:
                            print(name, str(length), 'errOutput: ', erroutput.strip('\n'), sep=', ', end='***\n')
                            print(name, str(length), 'stdOutput: ', stdoutput.strip('\n'), sep=', ', end='***\n')
                            if process.returncode < 0:
                                print('The process calling the aaRepeatScan program was killed by the system!')
                            print('chenhuilong1\n')
                            failFile.write(name + '\t' + str(length) + '\n')
                            print('***', file=failLogFile)
                            print(process.returncode, name, str(length), sep='\t', file=failLogFile)
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
                                    invertedPositionBigLengthLogFile.write(name + '\n')

                                invertedPositionFile = open(customFolderPath + '/' + name + '_positions_inverted.txt',
                                                            'w',
                                                            encoding='utf-8')
                                invertedPositionFileOriginal = open(
                                    customFolderPathOriginal + '/' + name + '_positions_inverted.txt', 'w',
                                    encoding='utf-8')

                                for chl in posTable.values:
                                    if chl[1] > chl[0] and chl[2] >= repeatLen:
                                        invertedPositionFile.write(
                                            str(int(chl[0])) + '\t' + str(int(chl[0] + chl[2] - 1)) + '\t' + str(int(chl[1])) + '\t'
                                            + str(int(chl[1] - chl[2] + 1)) + '\t' + str(int(chl[2])) + '\t' + str(
                                                '{:6f}'.format(chl[3]))
                                            + '\t' + str(int(chl[4])) + '\t' + str('{:6f}'.format(chl[5])) + '\t' + str(
                                                int(chl[6])) + '\n')
                                    if chl[2] >= repeatLen:
                                        invertedPositionFileOriginal.write(
                                            str(int(chl[0])) + '\t' + str(int(chl[0] + chl[2] - 1)) + '\t' + str(int(chl[1])) + '\t'
                                            + str(int(chl[1] - chl[2] + 1)) + '\t' + str(int(chl[2])) + '\t' + str(
                                                '{:6f}'.format(chl[3]))
                                            + '\t' + str(int(chl[4])) + '\t' + str('{:6f}'.format(chl[5])) + '\t' + str(
                                                int(chl[6])) + '\n')
                                invertedPositionFile.close()
                                invertedPositionFileOriginal.close()
                            except BaseException as e:
                                print(name, str(length), e, e.__traceback__.tb_lineno, sep='***')
                                print('chenhuilong2\n')
                                failFile.write(name + '\t' + str(length) + '\n')
                                print(name, str(length), e, e.__traceback__.tb_lineno, sep='***', file=failLogFile)

                                if os.path.exists(customFolderPath + '/' + name + '_positions_inverted.txt'):
                                    os.remove(customFolderPath + '/' + name + '_positions_inverted.txt')
                                if os.path.exists(customFolderPathOriginal + '/' + name + '_positions_inverted.txt'):
                                    os.remove(customFolderPathOriginal + '/' + name + '_positions_inverted.txt')
    
                if '2' in modeList:
                    if 'DNA' == fileType:
                        tempSingleInputFile1 = open('./temp.single.input.fasta1.txt', 'w', encoding='utf-8')
                        tempSingleInputFile1.write(seq)
                        tempSingleInputFile1.close()
                        tempSingleInputFile2 = open('./temp.single.input.fasta2.txt', 'w', encoding='utf-8')
                        tempSingleInputFile2.write(seq)
                        tempSingleInputFile2.close()
                        mode = 2
                        exactMatch = -1

                        command = "./bpRepeatScan {} {} {} {}".format(str(identityThr), str(mode), str(DNA_Matrix), str(exactMatch))
                        process = subprocess.Popen(command, universal_newlines=True, stdin=subprocess.PIPE,
                                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,
                                                   preexec_fn=os.setsid)
                        stdoutput, erroutput = process.communicate()

                        print("mode 2, {}, process.returncode: {}".format(name, process.returncode))
                        if process.returncode:
                            print(name, str(length), 'errOutput: ', erroutput.strip('\n'), sep=', ', end='***\n')
                            print(name, str(length), 'stdOutput: ', stdoutput.strip('\n'), sep=', ', end='***\n')
                            if process.returncode < 0:
                                print('The process calling the bpRepeatScan program was killed by the system!')
                            print('chenhuilong1\n')
                            failFile.write(name + '\t' + str(length) + '\n')
                            print('***', file=failLogFile)
                            print(process.returncode, name, str(length), sep='\t', file=failLogFile)
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
                                    reverseComplementPositionBigLengthLogFile.write(name + '\n')

                                reverseComplementPositionFile = open(
                                    customFolderPath + '/' + name + '_positions_reverse_complement.txt', 'w',
                                    encoding='utf-8')
                                reverseComplementPositionFileOriginal = open(
                                    customFolderPathOriginal + '/' + name + '_positions_reverse_complement.txt', 'w',
                                    encoding='utf-8')

                                for chl in posTable.values:
                                    if chl[1] > chl[0] and chl[2] >= repeatLen:
                                        reverseComplementPositionFile.write(
                                            str(int(chl[0])) + '\t' + str(int(chl[0] + chl[2] - 1)) + '\t' + str(int(chl[1])) + '\t'
                                            + str(int(chl[1] - chl[2] + 1)) + '\t' + str(int(chl[2])) + '\t' + str(
                                                '{:6f}'.format(chl[3]))
                                            + '\t' + str(int(chl[4])) + '\t' + str(int(chl[5])) + '\n')
                                    if chl[2] >= repeatLen:
                                        reverseComplementPositionFileOriginal.write(
                                            str(int(chl[0])) + '\t' + str(int(chl[0] + chl[2] - 1)) + '\t' + str(int(chl[1])) + '\t'
                                            + str(int(chl[1] - chl[2] + 1)) + '\t' + str(int(chl[2])) + '\t' + str(
                                                '{:6f}'.format(chl[3]))
                                            + '\t' + str(int(chl[4])) + '\t' + str(int(chl[5])) + '\n')
                                reverseComplementPositionFile.close()
                                reverseComplementPositionFileOriginal.close()
                            except BaseException as e:
                                print(name, str(length), e, e.__traceback__.tb_lineno, sep='***')
                                print('chenhuilong2\n')
                                failFile.write(name + '\t' + str(length) + '\n')
                                print(name, str(length), e, e.__traceback__.tb_lineno, sep='***', file=failLogFile)

                                if os.path.exists(customFolderPath + '/' + name + '_positions_reverse_complement.txt'):
                                    os.remove(customFolderPath + '/' + name + '_positions_reverse_complement.txt')
                                if os.path.exists(customFolderPathOriginal + '/' + name + '_positions_reverse_complement.txt'):
                                    os.remove(customFolderPathOriginal + '/' + name + '_positions_reverse_complement.txt')
            except BaseException as e:
                print(name, str(length), e, e.__traceback__.tb_lineno, sep='***')
                print('chenhuilong3\n')
                failFile.write(name + '\t' + str(length) + '\n')
                print(name, str(length), e, e.__traceback__.tb_lineno, sep='***', file=failLogFile)

            removeThreeFiles()

        logFile.close()
        failFile.close()
        failLogFile.close()
        if '0' in modeList:
            directPositionBigLengthLogFile.close()
        if '1' in modeList:
            invertedPositionBigLengthLogFile.close()
        if '2' in modeList:
            reverseComplementPositionBigLengthLogFile.close()
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
        batchExportPosition(chlFasta,exportFolderPositionPath,fileType,DNA_Matrix,aminoAcidMatrix,modeList,identityThr,similarityThr,repeatLen)

        totalNumber = len(list(chlFasta.keys()))
        
        folder = '/positions'
        customFailFilePath = exportFolderPositionPath + folder + '/failure_sequences.txt'
        with open(customFailFilePath, 'r', encoding='utf-8') as failFile:
            deduplicationFailNameDict = {}
            for line in failFile:
                print(line.strip('\n'))
                lineList = line.strip('\n').split('\t')
                deduplicationFailNameDict[lineList[0]] = ''
            failNameList = list(deduplicationFailNameDict.keys())
            failNumber = len(failNameList)

        # 'positions' folder
        customFolderPath = exportFolderPositionPath + folder
        customFileList = os.listdir(customFolderPath)

        deduplicationDict = {}
        for fileName in customFileList:
            if '_positions_' in fileName:
                name = re.findall("(.*?)_positions_", fileName)[0]
                deduplicationDict[name] = ''
                if name in failNameList:
                    print("For 'positions' folder, This '{}' is in 'failure_sequences.txt'. it's illogical!".format(name))
                    os.remove(customFolderPath + '/' + fileName)

        successNumber = len(list(deduplicationDict.keys()))

        customFileList = os.listdir(customFolderPath)
        deduplicationAfterDeleteDict = {}
        for fileName in customFileList:
            if '_positions_' in fileName:
                name = re.findall("(.*?)_positions_", fileName)[0]
                deduplicationAfterDeleteDict[name] = ''

        successNumberAfterDelete = len(list(deduplicationAfterDeleteDict.keys()))

        print("For 'positions',")
        print("This fasta file with {} sequences has a total of {} failures and {} successes this time!".format(
            totalNumber, failNumber, successNumber))
        print(
            "If a successful ID appears in 'failure_sequences.txt', the number result after deleting the position file corresponding to this ID is: \n "
            "This fasta file with {} sequences has a total of {} failures and {} successes this time!".format(
                totalNumber, failNumber, successNumberAfterDelete))

        # positions_original
        folderOriginal = '/positions_original'
        customFolderPathOriginal = exportFolderPositionPath + folderOriginal
        customOriginalFileList = os.listdir(customFolderPathOriginal)
        while 'all_sequences_length.txt' in customOriginalFileList:
            customOriginalFileList.remove('all_sequences_length.txt')

        deduplicationDictOriginal = {}
        for fileNameOriginal in customOriginalFileList:
            name = re.findall("(.*?)_positions_", fileNameOriginal)[0]
            deduplicationDictOriginal[name] = ''
            if name in failNameList:
                print("For 'positions_original' folder, This '{}' is in 'failure_sequences.txt'. it's illogical!".format(name))
                os.remove(customFolderPathOriginal + '/' + fileNameOriginal)

        successNumberOriginal = len(list(deduplicationDictOriginal.keys()))

        customOriginalFileList = os.listdir(customFolderPathOriginal)
        while 'all_sequences_length.txt' in customOriginalFileList:
            customOriginalFileList.remove('all_sequences_length.txt')

        deduplicationAfterDeleteDictOriginal = {}
        for fileNameOriginal in customOriginalFileList:
            name = re.findall("(.*?)_positions_", fileNameOriginal)[0]
            deduplicationAfterDeleteDictOriginal[name] = ''

        successNumberAfterDeleteOriginal = len(list(deduplicationAfterDeleteDictOriginal.keys()))

        print("For 'positions_original',")
        print("This fasta file with {} sequences has a total of {} failures and {} successes this time!".format(totalNumber, failNumber, successNumberOriginal))
        print("If a successful ID appears in 'failure_sequences.txt', the number result after deleting the position file corresponding to this ID is: \n "
              "This fasta file with {} sequences has a total of {} failures and {} successes this time!".format(totalNumber, failNumber, successNumberAfterDeleteOriginal))

        if totalNumber == failNumber + successNumber and totalNumber == failNumber + successNumberOriginal:
            print("For 'positions' and 'positions_original', the sum of this number of successes and failures is correct!")
        else:
            state = -1

        # 'positions' folder
        intersectionList = list(set(failNameList) & set(list(deduplicationDict.keys())))
        unionList = list(set(failNameList) | set(list(deduplicationDict.keys())))
        besidesIntersectionList = list(set(unionList) ^ set(list(chlFasta.keys())))
            
        # positions_original
        intersectionListOriginal = list(set(failNameList) & set(list(deduplicationDictOriginal.keys())))
        unionListOriginal = list(set(failNameList) | set(list(deduplicationDictOriginal.keys())))
        besidesIntersectionListOriginal = list(set(unionListOriginal) ^ set(list(chlFasta.keys())))

        if [] == intersectionList and [] == besidesIntersectionList and [] == intersectionListOriginal and [] == besidesIntersectionListOriginal and totalNumber == failNumber + successNumber and totalNumber == failNumber + successNumberOriginal:
            print(
                "For 'positions' and 'positions_original', the sum of this number of successes and failures and the ID is correct!\n")
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

# created by Huilong Chen, May 21, 2022!
# revised by Huilong Chen, May 27, 2022!
# revised by Huilong Chen, July 23, 2022! A new CyDotian.
# revised by Huilong Chen, July 25, 2022! Complete the code of other modes.
# revised by Huilong Chen, July 26, 2022! Optimize.
# revised by Huilong Chen, August 1, 2022! Optimize.
