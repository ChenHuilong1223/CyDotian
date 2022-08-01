#define _CRT_SECURE_NO_WARNINGS 1
#pragma warning(disable:6386)
#pragma warning(disable:6385)
#pragma warning(disable:4267)
#pragma warning(disable:4477)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "bpRepeatScan.h"
void ReverseComplement(char* dest, char* sour, int strlen)
{
	int i;
	for (i = strlen - 1;i >= 0;i--)
	{
		switch (sour[i])
		{
		case 'A':
			dest[strlen - 1 - i] = 'T';
			break;
		case 'T':
			dest[strlen - 1 - i] = 'A';
			break;
		case 'C':
			dest[strlen - 1 - i] = 'G';
			break;
		case 'G':
			dest[strlen - 1 - i] = 'C';
			break;
		case 'N':
			dest[strlen - 1 - i] = 'N';
			break;
		default:
			break;
		}
	}
	dest[strlen - 1 - i] = '\0';
}
void Reverse(char* dest, char* sour, int strlen)
{
	int i;
	for (i = 0; i < strlen; i++)
	{
		dest[i] = sour[strlen - 1 - i];
	}
	dest[i] = '\0';
}
int Cmp_scoreSite_by_score(const void* p1, const void* p2)
{
	return ((scoreSite_*)p2)->score - ((scoreSite_*)p1)->score;
}
int Cmp_scoreSite_by_site0(const void* p1, const void* p2)
{
	return ((scoreSite_*)p2)->site[0] - ((scoreSite_*)p1)->site[0];
}

int ExactMatch(char a, char b)
{
	if (a == b)
	{
		return 5;
	}
	else
	{
		return -4;
	}
}
int FindMatrixIndex(char c)
{
	switch (c)
	{
	case 'A':
		return 0;
		break;
	case 'T':
		return 1;
		break;
	case 'C':
		return 2;
		break;
	case 'G':
		return 3;
		break;
	case 'N':
		return 4;
		break;
	default:
		return -1;
		break;
	}
}
int NonExactMatch(char a, char b, int matrix)
{
	int i, j;
	switch (matrix)
	{
	case 0:
	{
		int BLASTMatrix[5][5] = {
			5,-4,-4,-4,-2,
			-4,5,-4,-4,-2,
			-4,-4,5,-4,-2,
			-4,-4,-4,5,-2,
			-2,-2,-2,-2,-1};
		i = FindMatrixIndex(a);
		j = FindMatrixIndex(b);
		return BLASTMatrix[i][j];
		break;
	}
	case 1:
	{
		int transitionTransversionMatrix[5][5] = {
			1,-5,-5,-1,-3,
			-5,1,-1,-5,-3,
			-5,-1,1,-5,-3,
			-1,-5,-5,1,-3,
			-3,-3,-3,-3,-2};
		i = FindMatrixIndex(a);
		j = FindMatrixIndex(b);
		return transitionTransversionMatrix[i][j];
		break;
	}
	default:
		return 0;
		break;
	}
}
strace_ AlgoScore(int diagonal, int match)
{
	if (diagonal + match > 0)
	{
		strace_ fsy = { 'F', diagonal + match };
		return fsy;
	}
	else
	{
		strace_ fsy = { '\0',0 };
		return fsy;
	}
}

void ReverseComplementFun(char* seqArr1, char* seqArr2, double ideThr, int matrix, int exact)
{
	int row, col;
	int i = 0, j = 0;
	int n, a = 0;
	int k = 1, o = 1, l = 0;
	int count;
	int diago, match;
	int getScore;

	double identity;
	char* revCompSeqArr;
	int** scores = NULL;
	char** directs = NULL;
	scoreSite_* scoPos;
	strace_ chl = { '\0',0 };

	row = strlen(seqArr1) + 1;
	col = strlen(seqArr2) + 1;
	int scoreSiteSize = (row - 1) * (col - 1);

	if (row <= 0 || col <= 0)
	{
		printf("The number of rows or columns is illegal!\n");
		exit(-1);
	}

	scores = (int**)malloc(row * sizeof(int*));
	if (NULL == scores) {
		printf("Unable to dynamically apply memory!\n");
		exit(-1);
	}
	for (int i = 0; i < row;i++)
	{
		scores[i] = (int*)malloc(col * sizeof(int));
		if (NULL == scores[i])
		{
			printf("Unable to dynamically apply memory!\n");
			exit(-1);
		}
	}

	directs = (char**)malloc(row * sizeof(char*));
	if (NULL == directs) {
		printf("Unable to dynamically apply memory!\n");
		exit(-1);
	}
	for (int i = 0; i < row;i++)
	{
		directs[i] = (char*)malloc(col * sizeof(char));
		if (NULL == directs[i])
		{
			printf("Unable to dynamically apply memory!\n");
			exit(-1);
		}
	}

	scoPos = (scoreSite_*)malloc(scoreSiteSize * sizeof(scoreSite_));
	if (NULL == scoPos) {
		printf("Unable to dynamically apply memory!\n");
		exit(-1);
	}

	revCompSeqArr = (char*)malloc(col * sizeof(char));
	if (NULL == revCompSeqArr)
	{
		printf("Unable to dynamically apply memory!\n");
		exit(-1);
	}
	ReverseComplement(revCompSeqArr, seqArr2, col - 1);

	for (i = 0;i < row;i++)
	{
		scores[i][0] = 0;
		directs[i][0] = '\0';
	}

	for (j = 0;j < col;j++)
	{
		scores[0][j] = 0;
		directs[0][j] = '\0';
	}

	n = 0;
	if (0 == exact)
	{
		for (i = 1;i < row;i++)
		{
			for (j = 1;j < col;j++)
			{
				diago = scores[i - 1][j - 1];
				match = ExactMatch(seqArr1[i - 1], revCompSeqArr[j - 1]);
				chl = AlgoScore(diago, match);
				scores[i][j] = chl.score;
				directs[i][j] = chl.direction;
				scoPos[n].score = chl.score;
				scoPos[n].site[0] = i;
				scoPos[n].site[1] = j;
				n++;
			}
		}
	}
	else
	{
		for (i = 1;i < row;i++)
		{
			for (j = 1;j < col;j++)
			{
				diago = scores[i - 1][j - 1];
				match = NonExactMatch(seqArr1[i - 1], revCompSeqArr[j - 1], matrix);
				chl = AlgoScore(diago, match);
				scores[i][j] = chl.score;
				directs[i][j] = chl.direction;
				scoPos[n].score = chl.score;
				scoPos[n].site[0] = i;
				scoPos[n].site[1] = j;
				n++;
			}
		}
	}

	qsort(scoPos, scoreSiteSize, sizeof(scoPos[0]), Cmp_scoreSite_by_score);

	int state = -1;
	int num = 1;
	for (n = 0; n < scoreSiteSize; n++)
	{
		if (0 == scoPos[n].score)
		{
			break;
		}
		if (scoPos[n].score != state)
		{
			state = scoPos[n].score;
			if (num != 1)
			{
				qsort(&scoPos[n - num], num, sizeof(scoPos[0]), Cmp_scoreSite_by_site0);
			}
			num = 1;
		}
		else
		{
			num++;
		}
	}
	if (num != 1)
	{
		qsort(&scoPos[n - num], num, sizeof(scoPos[0]), Cmp_scoreSite_by_site0);
	}

	FILE* pfWrite = fopen("position.txt", "w");
	count = 0;
	getScore = 0;
	if (0 == exact)
	{
		for (n = 0;n < scoreSiteSize;n++)
		{
			i = scoPos[n].site[0];
			j = scoPos[n].site[1];
			if (0 == scores[i][j])
			{
				;
			}
			else
			{
				if (seqArr1[i - 1] == revCompSeqArr[j - 1])
				{
					while (k)
					{
						if ('F' == directs[i][j] && 5 == scores[i][j] - scores[i - 1][j - 1])
						{
							scores[i][j] = 0;
							directs[i][j] = '\0';
							i--;
							j--;
							l++;
						}
						else
						{
							if (seqArr1[i + 1 - 1] == revCompSeqArr[j + 1 - 1])
							{
								if (l >= 1)
								{
									fprintf(pfWrite, "%d\t%d\t%d\n", i+1, col-1-j, l);
								}
							}
							else
							{
								while (o)
								{
									i++;
									j++;
									l--;
									if (seqArr1[i + 1 - 1] == revCompSeqArr[j + 1 - 1])
									{
										if (l >= 1)
										{
											fprintf(pfWrite, "%d\t%d\t%d\n", i+1, col-1-j, l);
										}
										break;
									}
								}
							}
							count = 0;
							l = 0;
							break;
						}
					}
				}
			}
		}
	}
	else
	{
		for (n = 0;n < scoreSiteSize;n++)
		{
			i = scoPos[n].site[0];
			j = scoPos[n].site[1];
			if (0 == scores[i][j])
			{
				;
			}
			else
			{
				if (seqArr1[i - 1] == revCompSeqArr[j - 1])
				{
					while (k)
					{
						if ('F' == directs[i][j])
						{
							if (seqArr1[i - 1] == revCompSeqArr[j - 1])
							{
								count++;
							}
							getScore = getScore + NonExactMatch(seqArr1[i - 1], revCompSeqArr[j - 1], matrix);
							scores[i][j] = 0;
							directs[i][j] = '\0';
							i--;
							j--;
							l++;
						}
						else
						{
							if (seqArr1[i + 1 - 1] == revCompSeqArr[j + 1 - 1])
							{
								identity = (double)count / (double)l;
								if (identity >= ideThr && l >= 1)
								{
									fprintf(pfWrite, "%d\t%d\t%d\t%lf\t%d\t%d\n", i+1, col-1-j, l, identity, l-count, getScore);
								}
							}
							else
							{
								while (o)
								{
									getScore = getScore - NonExactMatch(seqArr1[i + 1 - 1], revCompSeqArr[j + 1 - 1], matrix);
									i++;
									j++;
									l--;
									if (seqArr1[i + 1 - 1] == revCompSeqArr[j + 1 - 1])
									{
										identity = (double)count / (double)l;
										if (identity >= ideThr && l >= 1)
										{
											fprintf(pfWrite, "%d\t%d\t%d\t%lf\t%d\t%d\n", i+1, col-1-j, l, identity, l-count, getScore);
										}
										break;
									}
								}
							}
							count = 0;
							l = 0;
							getScore = 0;
							break;
						}
					}
				}
			}
		}
	}

	fclose(pfWrite);

	for (int i = 0;i < row; i++)
	{
		free(scores[i]);
		free(directs[i]);
	}

	free(scores);
	free(directs);
	free(scoPos);
	free(revCompSeqArr);
}

void ReverseFun(char* seqArr1, char* seqArr2, double ideThr, int matrix, int exact)
{
	int row, col;
	int i = 0, j = 0;
	int n, a = 0;
	int k = 1, o = 1, l = 0;
	int count;
	int diago, match;
	int getScore;

	double identity;
	char* revSeqArr;
	int** scores = NULL;
	char** directs = NULL;
	scoreSite_* scoPos;
	strace_ chl = { '\0',0 };

	row = strlen(seqArr1) + 1;
	col = strlen(seqArr2) + 1;
	int scoreSiteSize = (row - 1) * (col - 1);

	if (row <= 0 || col <= 0)
	{
		printf("The number of rows or columns is illegal!\n");
		exit(-1);
	}

	scores = (int**)malloc(row * sizeof(int*));
	if (NULL == scores) {
		printf("Unable to dynamically apply memory!\n");
		exit(-1);
	}
	for (int i = 0; i < row;i++)
	{
		scores[i] = (int*)malloc(col * sizeof(int));
		if (NULL == scores[i])
		{
			printf("Unable to dynamically apply memory!\n");
			exit(-1);
		}
	}

	directs = (char**)malloc(row * sizeof(char*));
	if (NULL == directs) {
		printf("Unable to dynamically apply memory!\n");
		exit(-1);
	}
	for (int i = 0; i < row;i++)
	{
		directs[i] = (char*)malloc(col * sizeof(char));
		if (NULL == directs[i])
		{
			printf("Unable to dynamically apply memory!\n");
			exit(-1);
		}
	}
	scoPos = (scoreSite_*)malloc(scoreSiteSize * sizeof(scoreSite_));
	if (NULL == scoPos) {
		printf("Unable to dynamically apply memory!\n");
		exit(-1);
	}
	revSeqArr = (char*)malloc(col * sizeof(char));
	if (NULL == revSeqArr)
	{
		printf("Unable to dynamically apply memory!\n");
		exit(-1);
	}
	Reverse(revSeqArr, seqArr2, col - 1);

	for (i = 0;i < row;i++)
	{
		scores[i][0] = 0;
		directs[i][0] = '\0';
	}

	for (j = 0;j < col;j++)
	{
		scores[0][j] = 0;
		directs[0][j] = '\0';
	}

	n = 0;
	if (0 == exact)
	{
		for (i = 1;i < row;i++)
		{
			for (j = 1;j < col;j++)
			{
				diago = scores[i - 1][j - 1];
				match = ExactMatch(seqArr1[i - 1], revSeqArr[j - 1]);
				chl = AlgoScore(diago, match);
				scores[i][j] = chl.score;
				directs[i][j] = chl.direction;
				scoPos[n].score = chl.score;
				scoPos[n].site[0] = i;
				scoPos[n].site[1] = j;
				n++;
			}
		}
	}
	else
	{
		for (i = 1;i < row;i++)
		{
			for (j = 1;j < col;j++)
			{
				diago = scores[i - 1][j - 1];
				match = NonExactMatch(seqArr1[i - 1], revSeqArr[j - 1], matrix);
				chl = AlgoScore(diago, match);
				scores[i][j] = chl.score;
				directs[i][j] = chl.direction;
				scoPos[n].score = chl.score;
				scoPos[n].site[0] = i;
				scoPos[n].site[1] = j;
				n++;
			}
		}
	}

	qsort(scoPos, scoreSiteSize, sizeof(scoPos[0]), Cmp_scoreSite_by_score);

	int state = -1;
	int num = 1;
	for (n = 0; n < scoreSiteSize; n++)
	{
		if (0 == scoPos[n].score)
		{
			break;
		}
		if (scoPos[n].score != state)
		{
			state = scoPos[n].score;
			if (num != 1)
			{
				qsort(&scoPos[n - num], num, sizeof(scoPos[0]), Cmp_scoreSite_by_site0);
			}
			num = 1;
		}
		else
		{
			num++;
		}
	}
	if (num != 1)
	{
		qsort(&scoPos[n - num], num, sizeof(scoPos[0]), Cmp_scoreSite_by_site0);
	}

	FILE* pfWrite = fopen("position.txt", "w");
	count = 0;
	getScore = 0;
	if (0 == exact)
	{
		for (n = 0;n < scoreSiteSize;n++)
		{
			i = scoPos[n].site[0];
			j = scoPos[n].site[1];
			if (0 == scores[i][j])
			{
				;
			}
			else
			{
				if (seqArr1[i - 1] == revSeqArr[j - 1])
				{
					while (k)
					{
						if ('F' == directs[i][j] && 5 == scores[i][j] - scores[i - 1][j - 1])
						{
							scores[i][j] = 0;
							directs[i][j] = '\0';
							i--;
							j--;
							l++;
						}
						else
						{
							if (seqArr1[i + 1 - 1] == revSeqArr[j + 1 - 1])
							{
								if (l >= 1)
								{
									fprintf(pfWrite, "%d\t%d\t%d\n", i+1, col-1-j, l);
								}
							}
							else
							{
								while (o)
								{
									i++;
									j++;
									l--;
									if (seqArr1[i + 1 - 1] == revSeqArr[j + 1 - 1])
									{
										if (l >= 1)
										{
											fprintf(pfWrite, "%d\t%d\t%d\n", i+1, col-1-j, l);
										}
										break;
									}
								}
							}
							count = 0;
							l = 0;
							break;
						}
					}
				}
			}
		}
	}
	else
	{
		for (n = 0;n < scoreSiteSize;n++)
		{
			i = scoPos[n].site[0];
			j = scoPos[n].site[1];
			if (0 == scores[i][j])
			{
				;
			}
			else
			{
				if (seqArr1[i - 1] == revSeqArr[j - 1])
				{
					while (k)
					{
						if ('F' == directs[i][j])
						{
							if (seqArr1[i - 1] == revSeqArr[j - 1])
							{
								count++;
							}
							getScore = getScore + NonExactMatch(seqArr1[i - 1], revSeqArr[j - 1], matrix);
							scores[i][j] = 0;
							directs[i][j] = '\0';
							i--;
							j--;
							l++;
						}
						else
						{
							if (seqArr1[i + 1 - 1] == revSeqArr[j + 1 - 1])
							{
								identity = (double)count / (double)l;
								if (identity >= ideThr && l >= 1)
								{
									fprintf(pfWrite, "%d\t%d\t%d\t%lf\t%d\t%d\n", i+1, col-1-j, l, identity, l-count, getScore);
								}
							}
							else
							{
								while (o)
								{
									getScore = getScore - NonExactMatch(seqArr1[i + 1 - 1], revSeqArr[j + 1 - 1], matrix);
									i++;
									j++;
									l--;
									if (seqArr1[i + 1 - 1] == revSeqArr[j + 1 - 1])
									{
										identity = (double)count / (double)l;
										if (identity >= ideThr && l >= 1)
										{
											fprintf(pfWrite, "%d\t%d\t%d\t%lf\t%d\t%d\n", i+1, col-1-j, l, identity, l-count, getScore);
										}
										break;
									}
								}
							}
							count = 0;
							l = 0;
							getScore = 0;
							break;
						}
					}
				}
			}
		}
	}

	fclose(pfWrite);

	for (int i = 0;i < row; i++)
	{
		free(scores[i]);
		free(directs[i]);
	}
	free(scores);
	free(directs);
	free(scoPos);
	free(revSeqArr);
}

void DirectFun(char* seqArr1, char* seqArr2, double ideThr, int matrix, int exact)
{

	int row, col;
	int i = 0, j = 0;
	int n, a = 0;
	int k = 1, o = 1, l = 0;
	int count;
	int diago, match;
	int getScore;
	double identity;

	int** scores = NULL;
	char** directs = NULL;
	scoreSite_* scoPos;
	strace_ chl = { '\0',0 };

	row = strlen(seqArr1) + 1;
	col = strlen(seqArr2) + 1;
	int scoreSiteSize = (row - 1) * (col - 1);

	if (row <= 0 || col <= 0)
	{
		printf("The number of rows or columns is illegal!\n");
		exit(-1);
	}

	scores = (int**)malloc(row * sizeof(int*));
	if (NULL == scores) {
		printf("Unable to dynamically apply memory!\n");
		exit(-1);
	}
	for (int i = 0; i < row;i++)
	{
		scores[i] = (int*)malloc(col * sizeof(int));
		if (NULL == scores[i])
		{
			printf("Unable to dynamically apply memory!\n");
			exit(-1);
		}
	}
	directs = (char**)malloc(row * sizeof(char*));
	if (NULL == directs) {
		printf("Unable to dynamically apply memory!\n");
		exit(-1);
	}
	for (int i = 0; i < row;i++)
	{
		directs[i] = (char*)malloc(col * sizeof(char));
		if (NULL == directs[i])
		{
			printf("Unable to dynamically apply memory!\n");
			exit(-1);
		}
	}
	scoPos = (scoreSite_*)malloc(scoreSiteSize * sizeof(scoreSite_));
	if (NULL == scoPos) {
		printf("Unable to dynamically apply memory!\n");
		exit(-1);
	}

	for (i = 0;i < row;i++)
	{
		scores[i][0] = 0;
		directs[i][0] = '\0';
	}

	for (j = 0;j < col;j++)
	{
		scores[0][j] = 0;
		directs[0][j] = '\0';
	}

	n = 0;
	if (0 == exact)
	{
		for (i = 1;i < row;i++)
		{
			for (j = 1;j < col;j++)
			{
				diago = scores[i - 1][j - 1];
				match = ExactMatch(seqArr1[i - 1], seqArr2[j - 1]);
				chl = AlgoScore(diago, match);
				scores[i][j] = chl.score;
				directs[i][j] = chl.direction;
				scoPos[n].score = chl.score;
				scoPos[n].site[0] = i;
				scoPos[n].site[1] = j;
				n++;
			}
		}
	}
	else
	{
		for (i = 1;i < row;i++)
		{
			for (j = 1;j < col;j++)
			{
				diago = scores[i - 1][j - 1];
				match = NonExactMatch(seqArr1[i - 1], seqArr2[j - 1], matrix);
				chl = AlgoScore(diago, match);
				scores[i][j] = chl.score;
				directs[i][j] = chl.direction;
				scoPos[n].score = chl.score;
				scoPos[n].site[0] = i;
				scoPos[n].site[1] = j;
				n++;
			}
		}
	}
	
	qsort(scoPos, scoreSiteSize, sizeof(scoPos[0]), Cmp_scoreSite_by_score);

	int state = -1;
	int num = 1;
	for (n = 0; n < scoreSiteSize; n++)
	{
		if (0 == scoPos[n].score)
		{
			break;
		}
		if (scoPos[n].score != state)
		{
			state = scoPos[n].score;
			if (num != 1)
			{
				qsort(&scoPos[n - num], num, sizeof(scoPos[0]), Cmp_scoreSite_by_site0);
			}
			num = 1;
		}
		else
		{
			num++;
		}
	}
	if (num != 1)
	{
		qsort(&scoPos[n - num], num, sizeof(scoPos[0]), Cmp_scoreSite_by_site0);
	}

	FILE* pfWrite = fopen("position.txt", "w");
	count = 0;
	getScore = 0;
	if (0 == exact)
	{
		for (n = 0;n < scoreSiteSize;n++)
		{
			i = scoPos[n].site[0];
			j = scoPos[n].site[1];
			if (0 == scores[i][j])
			{
				;
			}
			else
			{
				if (seqArr1[i - 1] == seqArr2[j - 1])
				{
					while (k)
					{
						if ('F' == directs[i][j] && 5 == scores[i][j] - scores[i - 1][j - 1])
						{
							scores[i][j] = 0;
							directs[i][j] = '\0';
							i--;
							j--;
							l++;
						}
						else
						{
							if (seqArr1[i + 1 - 1] == seqArr2[j + 1 - 1])
							{
								if (l >= 1)
								{
									fprintf(pfWrite, "%d\t%d\t%d\n", i+1, j+1, l);
								}
							}
							else
							{
								while (o)
								{
									i++;
									j++;
									l--;
									if (seqArr1[i + 1 - 1] == seqArr2[j + 1 - 1])
									{
										if (l >= 1)
										{
											fprintf(pfWrite, "%d\t%d\t%d\n", i+1, j+1, l);
										}
										break;
									}
								}
							}
							count = 0;
							l = 0;
							break;
						}
					}
				}
			}
		}
	}
	else
	{
		for (n = 0;n < scoreSiteSize;n++)
		{
			i = scoPos[n].site[0];
			j = scoPos[n].site[1];
			if (0 == scores[i][j])
			{
				;
			}
			else
			{
				if (seqArr1[i - 1] == seqArr2[j - 1])
				{
					while (k)
					{
						if ('F' == directs[i][j])
						{
							if (seqArr1[i - 1] == seqArr2[j - 1])
							{
								count++;
							}
							getScore = getScore + NonExactMatch(seqArr1[i - 1], seqArr2[j - 1], matrix);
							scores[i][j] = 0;
							directs[i][j] = '\0';
							i--;
							j--;
							l++;
						}
						else
						{
							if (seqArr1[i + 1 - 1] == seqArr2[j + 1 - 1])
							{
								identity = (double)count / (double)l;
								if (identity >= ideThr && l >= 1)
								{
									fprintf(pfWrite, "%d\t%d\t%d\t%lf\t%d\t%d\n", i+1, j+1, l, identity, l-count, getScore);
								}
							}
							else
							{
								while (o)
								{
									getScore = getScore - NonExactMatch(seqArr1[i + 1 - 1], seqArr2[j + 1 - 1], matrix);
									i++;
									j++;
									l--;
									if (seqArr1[i + 1 - 1] == seqArr2[j + 1 - 1])
									{
										identity = (double)count / (double)l;
										if (identity >= ideThr && l >= 1)
										{
											fprintf(pfWrite, "%d\t%d\t%d\t%lf\t%d\t%d\n", i+1, j+1, l, identity, l-count, getScore);
										}
										break;
									}
								}
							}
							count = 0;
							l = 0;
							getScore = 0;
							break;
						}
					}
				}
			}
		}
	}

	fclose(pfWrite);

	for (int i = 0;i < row; i++)
	{
		free(scores[i]);
		free(directs[i]);
	}
	free(scores);
	free(directs);
	free(scoPos);
}

int main(int argc, char* argv[])
{	
	int seqLengthIndex1 = 0;
	int seqLength1 = 0;
	int c1;
	FILE* fpSingleInputSeq1 = fopen("temp.single.input.fasta1.txt", "r");
	if(!fpSingleInputSeq1)
	{
		perror("File opening failed");
	}

	while ((c1 = fgetc(fpSingleInputSeq1)) != EOF)
	{
		seqLength1++;
	}
	if (ferror(fpSingleInputSeq1))
	{
		puts("I/O error when reading");
	}
	else if (feof(fpSingleInputSeq1))
	{
		;
	}

	char* seqArr1 = (char*)malloc((seqLength1+1) * sizeof(char));
	if (NULL == seqArr1)
	{
		printf("Unable to dynamically apply memory!\n");
		exit(-1);
	}

	rewind(fpSingleInputSeq1);

	while ((c1 = fgetc(fpSingleInputSeq1)) != EOF)
	{
		seqArr1[seqLengthIndex1++] = c1;
	}

	seqArr1[seqLengthIndex1] = '\0';

	if (ferror(fpSingleInputSeq1))
	{
		puts("I/O error when reading");
	}
	else if (feof(fpSingleInputSeq1))
	{
		;
	}

	fclose(fpSingleInputSeq1);

	int seqLengthIndex2 = 0;
	int seqLength2 = 0;
	int c2;
	FILE* fpSingleInputSeq2 = fopen("temp.single.input.fasta2.txt", "r");
	if(!fpSingleInputSeq2)
	{
		perror("File opening failed");
	}

	while ((c2 = fgetc(fpSingleInputSeq2)) != EOF)
	{
		seqLength2++;
	}
	if (ferror(fpSingleInputSeq2))
	{
		puts("I/O error when reading");
	}
	else if (feof(fpSingleInputSeq2))
	{
		;
	}	

	char* seqArr2 = (char*)malloc((seqLength2+1) * sizeof(char));
	if (NULL == seqArr2)
	{
		printf("Unable to dynamically apply memory!\n");
		exit(-1);
	}

	rewind(fpSingleInputSeq2);

	while ((c2 = fgetc(fpSingleInputSeq2)) != EOF)
	{
		seqArr2[seqLengthIndex2++] = c2;
	}

	seqArr2[seqLengthIndex2] = '\0';

	if (ferror(fpSingleInputSeq2))
	{
		puts("I/O error when reading");
	}
	else if (feof(fpSingleInputSeq2))
	{
		;
	}
	fclose(fpSingleInputSeq2);

	double ideThr = atof(argv[1]);
	int mode = atoi(argv[2]);
	int matrix = atoi(argv[3]);
	int exact = atoi(argv[4]);

	if (0 == mode)
	{
		DirectFun(seqArr1, seqArr2, ideThr, matrix, exact);
	}
	else if (1 == mode)
	{
		ReverseFun(seqArr1, seqArr2, ideThr, matrix, exact);
	}
	else if (2 == mode)
	{
		ReverseComplementFun(seqArr1, seqArr2, ideThr, matrix, exact);
	}


	free(seqArr1);
	free(seqArr2);
	
	return 0;

	// created by Huilong Chen!
}
