#ifndef _REPEATSCAN_H_
#define _REPEATSCAN_H_
int aaFindMatrixIndex(char c);
int aaNonExactMatch(char a, char b, int matrix);
double CalculateSimilarity(char* seq2, char* seq1, int windowSize, int matrix);
double CalculateIdentity(char* seq2, char* seq1, int windowSize);
void Reverse(char* dest, char* sour, int strlen);
void ReverseFun(char* seqArr1, char* seqArr2, double ideSimThr, int windowSize, int aaMatrix, int seqType);
void DirectFun(char* seqArr1, char* seqArr2, double ideSimThr, int windowSize, int aaMatrix, int seqType);
int main(int argc, char* argv[]);
#endif