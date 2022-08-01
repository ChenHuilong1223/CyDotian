#ifndef _REPEATSCAN_H_
#define _REPEATSCAN_H_
void Reverse(char* dest, char* sour, int strlen);
typedef struct scoreSite
{
	int score;
	int site[2];
} scoreSite_;
int Cmp_scoreSite_by_score(const void* p1, const void* p2);
int Cmp_scoreSite_by_site0(const void* p1, const void* p2);
typedef struct strace
{
	char direction;
	int score;
} strace_;
int ExactMatch(char a, char b);
int FindMatrixIndex(char c);
int NonExactMatch(char a, char b, int matrix);
strace_ AlgoScore(int diagonal, int match);
void ReverseFun(char* seqArr1, char* seqArr2, double simThr, int matrix, int exact);
void DirectFun(char* seqArr1, char* seqArr2, double simThr, int matrix, int exact);
int main(int argc, char* argv[]);
#endif