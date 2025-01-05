import os, time, sys, math, argparse
# sys.stdout = codecs.getwriter("utf-8")(sys.stdout.detach())

parser = argparse.ArgumentParser(description="check your sequences!")
parser.add_argument('-f', '--inputSequenceType', metavar='fileType', required=True, help='Input sequence type, 0 for DNA, 1 for amino acid')
parser.add_argument('-s', '--inputSequenceFolderPath', metavar='sequence', required=True, help='Input sequence folder path')
args = parser.parse_args()

startTime = time.time()
fileType = args.inputSequenceType
state = 0
try:
    importFolder = args.inputSequenceFolderPath
    importFileList = os.listdir(importFolder)
    importFolderPath = os.path.abspath(importFolder)

    def y(x):
        return '\n'.join([x[80*i:80*(i+1)] for i in range(math.ceil(len(x) / 80))])

    def replaceIllegalCharacter(seq):
        legalCharacterList = ['A','T','C','G','N']
        replacedSeq = ''
        for ch in seq:
            if ch not in legalCharacterList:
                ch = 'N'
            replacedSeq += ch

        return replacedSeq

    if '0' == fileType:
        allowedSet = set('ATCGN')
    elif '1' == fileType:
        allowedSet = set('ARNDCQEGHILKMFPSTWYVBZX*') # ABCDEFGHI KLMN PQRST VWXYZ

    for fileName1 in importFileList:
        userResponse = 'no'
        checkResult = 'normal'
        filePath1 = os.path.join(importFolderPath, fileName1)
        with open(filePath1, 'r', encoding='utf-8') as fastaFile:
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

        for name, seq in chlFasta.items():
            if set(seq) - allowedSet:
                print("The sequence of '{}' in file '{}' has illegal characters, please check!".format(name, fileName1))
                checkResult = 'issue'
                if '0' == fileType:
                    print("CyDotian provides a way to handle this by replacing all illegal characters with N, do you agree? yes/no")

                    reply = input("Please input your reply:")
                    while True:
                        if reply.strip().lower() == 'no':
                            userResponse = 'no'
                            break
                        elif reply.strip().lower() == 'yes':
                            userResponse = 'yes'
                            break
                        else:
                            reply = input("Please input yes or no:")

                break

        if '0' == fileType:
            if 'yes' == userResponse:
                importFolderParentPath = os.path.dirname(importFolderPath)
                replacedFolderPath = importFolderParentPath + '/' + os.path.basename(importFolderPath) + '_replace'

                isExists = os.path.exists(replacedFolderPath)
                if isExists:
                    pass
                else:
                    os.makedirs(replacedFolderPath, mode=0o777)

                replacedFilePath = replacedFolderPath + '/' + fileName1

                with open(replacedFilePath, 'w', encoding='utf-8') as replacedFile:
                    for name, seq in chlFasta.items():
                        replacedSeq = replaceIllegalCharacter(seq)
                        replacedFile.write(">{}\n{}\n".format(name, y(replacedSeq)))

                print("Your '{}' has been successfully replaced and stored in '{}'!".format(fileName1, replacedFilePath))
            elif 'no' == userResponse:
                pass

        if 'normal' == checkResult:
            print("Congratulations, your file '{}' can be analyzed by CyDotian!".format(fileName1))
except BaseException as e:
    state = -1
    print(e, e.__traceback__.tb_lineno)

if 0 == state:
    print('Please note that if you have too many N characters in your DNA sequence, a sequence that is inherently meaningless to analyze, please consider deleting it.'
          'Congratulations, the script worked and finished successfully!')
elif -1 == state:
    print('Sadly, the script did not complete properly, please check the output log, resolve the problem, and try again!')

endTime = time.time()
runTime = round(endTime - startTime)
hour = runTime//3600
minute = (runTime-3600*hour)//60
second = runTime-3600*hour-60*minute
print(f'The program running time: {hour}hour(s) {minute}minute(s) {second}second(s).')

# Created by Huilong Chen, July 15, 2022!
# revised by Huilong Chen, August 1, 2022! Optimize.
