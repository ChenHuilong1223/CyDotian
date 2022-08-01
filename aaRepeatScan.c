#define _CRT_SECURE_NO_WARNINGS 1
#pragma warning(disable:6386)
#pragma warning(disable:6385)
#pragma warning(disable:4267)
#pragma warning(disable:4477)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "aaRepeatScan.h"
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
	case 'R':
		return 1;
		break;
	case 'N':
		return 2;
		break;
	case 'D':
		return 3;
		break;
	case 'C':
		return 4;
		break;
	case 'Q':
		return 5;
		break;
	case 'E':
		return 6;
		break;
	case 'G':
		return 7;
		break;
	case 'H':
		return 8;
		break;
	case 'I':
		return 9;
		break;
	case 'L':
		return 10;
		break;
	case 'K':
		return 11;
		break;
	case 'M':
		return 12;
		break;
	case 'F':
		return 13;
		break;
	case 'P':
		return 14;
		break;
	case 'S':
		return 15;
		break;
	case 'T':
		return 16;
		break;
	case 'W':
		return 17;
		break;
	case 'Y':
		return 18;
		break;
	case 'V':
		return 19;
		break;
	case 'B':
		return 20;
		break;
	case 'Z':
		return 21;
		break;
	case 'X':
		return 22;
		break;
	case '*':
		return 23;
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
		int BLOSUM45[24][24] = {
		5,-2,-1,-2,-1,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-2,-2,0,-1,-1,-1,-5,
		-2,7,0,-1,-3,1,0,-2,0,-3,-2,3,-1,-2,-2,-1,-1,-2,-1,-2,-1,1,-1,-5,
		-1,0,6,2,-2,0,0,0,1,-2,-3,0,-2,-2,-2,1,0,-4,-2,-3,5,0,-1,-5,
		-2,-1,2,7,-3,0,2,-1,0,-4,-3,0,-3,-4,-1,0,-1,-4,-2,-3,6,1,-1,-5,
		-1,-3,-2,-3,12,-3,-3,-3,-3,-3,-2,-3,-2,-2,-4,-1,-1,-5,-3,-1,-2,-3,-1,-5,
		-1,1,0,0,-3,6,2,-2,1,-2,-2,1,0,-4,-1,0,-1,-2,-1,-3,0,4,-1,-5,
		-1,0,0,2,-3,2,6,-2,0,-3,-2,1,-2,-3,0,0,-1,-3,-2,-3,1,5,-1,-5,
		0,-2,0,-1,-3,-2,-2,7,-2,-4,-3,-2,-2,-3,-2,0,-2,-2,-3,-3,-1,-2,-1,-5,
		-2,0,1,0,-3,1,0,-2,10,-3,-2,-1,0,-2,-2,-1,-2,-3,2,-3,0,0,-1,-5,
		-1,-3,-2,-4,-3,-2,-3,-4,-3,5,2,-3,2,0,-2,-2,-1,-2,0,3,-3,-3,-1,-5,
		-1,-2,-3,-3,-2,-2,-2,-3,-2,2,5,-3,2,1,-3,-3,-1,-2,0,1,-3,-2,-1,-5,
		-1,3,0,0,-3,1,1,-2,-1,-3,-3,5,-1,-3,-1,-1,-1,-2,-1,-2,0,1,-1,-5,
		-1,-1,-2,-3,-2,0,-2,-2,0,2,2,-1,6,0,-2,-2,-1,-2,0,1,-2,-1,-1,-5,
		-2,-2,-2,-4,-2,-4,-3,-3,-2,0,1,-3,0,8,-3,-2,-1,1,3,0,-3,-3,-1,-5,
		-1,-2,-2,-1,-4,-1,0,-2,-2,-2,-3,-1,-2,-3,9,-1,-1,-3,-3,-3,-2,-1,-1,-5,
		1,-1,1,0,-1,0,0,0,-1,-2,-3,-1,-2,-2,-1,4,2,-4,-2,-1,0,0,-1,-5,
		0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-1,-1,2,5,-3,-1,0,0,-1,-1,-5,
		-2,-2,-4,-4,-5,-2,-3,-2,-3,-2,-2,-2,-2,1,-3,-4,-3,15,3,-3,-4,-2,-1,-5,
		-2,-1,-2,-2,-3,-1,-2,-3,2,0,0,-1,0,3,-3,-2,-1,3,8,-1,-2,-2,-1,-5,
		0,-2,-3,-3,-1,-3,-3,-3,-3,3,1,-2,1,0,-3,-1,0,-3,-1,5,-3,-3,-1,-5,
		-1,-1,5,6,-2,0,1,-1,0,-3,-3,0,-2,-3,-2,0,0,-4,-2,-3,5,1,-1,-5,
		-1,1,0,1,-3,4,5,-2,0,-3,-2,1,-1,-3,-1,0,-1,-2,-2,-3,1,5,-1,-5,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-5,
		-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,1
		};
		i = FindMatrixIndex(a);
		j = FindMatrixIndex(b);
		return BLOSUM45[i][j];
		break;
	}
	case 1:
	{
		int BLOSUM62[24][24] = {
		4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,-2,-1,0,-4,
		-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3,-1,0,-1,-4,
		-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3,3,0,-1,-4,
		-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3,4,1,-1,-4,
		0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4,
		-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2,0,3,-1,-4,
		-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4,
		0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3,-1,-2,-1,-4,
		-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3,0,0,-1,-4,
		-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3,-3,-3,-1,-4,
		-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1,-4,-3,-1,-4,
		-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2,0,1,-1,-4,
		-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1,-3,-1,-1,-4,
		-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1,-3,-3,-1,-4,
		-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2,-2,-1,-2,-4,
		1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2,0,0,0,-4,
		0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0,-1,-1,0,-4,
		-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3,-4,-3,-2,-4,
		-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1,-3,-2,-1,-4,
		0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4,-3,-2,-1,-4,
		-2,-1,3,4,-3,0,1,-1,0,-3,-4,0,-3,-3,-2,0,-1,-4,-3,-3,4,1,-1,-4,
		-1,0,0,1,-3,3,4,-2,0,-3,-3,1,-1,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4,
		0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,0,0,-2,-1,-1,-1,-1,-1,-4,
		-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,1,
		};
		i = FindMatrixIndex(a);
		j = FindMatrixIndex(b);
		return BLOSUM62[i][j];
		break;
	}
	case 2:
	{
		int BLOSUM80[24][24] = {
		5,-2,-2,-2,-1,-1,-1,0,-2,-2,-2,-1,-1,-3,-1,1,0,-3,-2,0,-2,-1,-1,-6,
		-2,6,-1,-2,-4,1,-1,-3,0,-3,-3,2,-2,-4,-2,-1,-1,-4,-3,-3,-1,0,-1,-6,
		-2,-1,6,1,-3,0,-1,-1,0,-4,-4,0,-3,-4,-3,0,0,-4,-3,-4,5,0,-1,-6,
		-2,-2,1,6,-4,-1,1,-2,-2,-4,-5,-1,-4,-4,-2,-1,-1,-6,-4,-4,5,1,-1,-6,
		-1,-4,-3,-4,9,-4,-5,-4,-4,-2,-2,-4,-2,-3,-4,-2,-1,-3,-3,-1,-4,-4,-1,-6,
		-1,1,0,-1,-4,6,2,-2,1,-3,-3,1,0,-4,-2,0,-1,-3,-2,-3,0,4,-1,-6,
		-1,-1,-1,1,-5,2,6,-3,0,-4,-4,1,-2,-4,-2,0,-1,-4,-3,-3,1,5,-1,-6,
		0,-3,-1,-2,-4,-2,-3,6,-3,-5,-4,-2,-4,-4,-3,-1,-2,-4,-4,-4,-1,-3,-1,-6,
		-2,0,0,-2,-4,1,0,-3,8,-4,-3,-1,-2,-2,-3,-1,-2,-3,2,-4,-1,0,-1,-6,
		-2,-3,-4,-4,-2,-3,-4,-5,-4,5,1,-3,1,-1,-4,-3,-1,-3,-2,3,-4,-4,-1,-6,
		-2,-3,-4,-5,-2,-3,-4,-4,-3,1,4,-3,2,0,-3,-3,-2,-2,-2,1,-4,-3,-1,-6,
		-1,2,0,-1,-4,1,1,-2,-1,-3,-3,5,-2,-4,-1,-1,-1,-4,-3,-3,-1,1,-1,-6,
		-1,-2,-3,-4,-2,0,-2,-4,-2,1,2,-2,6,0,-3,-2,-1,-2,-2,1,-3,-1,-1,-6,
		-3,-4,-4,-4,-3,-4,-4,-4,-2,-1,0,-4,0,6,-4,-3,-2,0,3,-1,-4,-4,-1,-6,
		-1,-2,-3,-2,-4,-2,-2,-3,-3,-4,-3,-1,-3,-4,8,-1,-2,-5,-4,-3,-2,-2,-1,-6,
		1,-1,0,-1,-2,0,0,-1,-1,-3,-3,-1,-2,-3,-1,5,1,-4,-2,-2,0,0,-1,-6,
		0,-1,0,-1,-1,-1,-1,-2,-2,-1,-2,-1,-1,-2,-2,1,5,-4,-2,0,-1,-1,-1,-6,
		-3,-4,-4,-6,-3,-3,-4,-4,-3,-3,-2,-4,-2,0,-5,-4,-4,11,2,-3,-5,-3,-1,-6,
		-2,-3,-3,-4,-3,-2,-3,-4,2,-2,-2,-3,-2,3,-4,-2,-2,2,7,-2,-3,-3,-1,-6,
		0,-3,-4,-4,-1,-3,-3,-4,-4,3,1,-3,1,-1,-3,-2,0,-3,-2,4,-4,-3,-1,-6,
		-2,-1,5,5,-4,0,1,-1,-1,-4,-4,-1,-3,-4,-2,0,-1,-5,-3,-4,5,0,-1,-6,
		-1,0,0,1,-4,4,5,-3,0,-4,-3,1,-1,-4,-2,0,-1,-3,-3,-3,0,5,-1,-6,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-6,
		-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,1
		};
		i = FindMatrixIndex(a);
		j = FindMatrixIndex(b);
		return BLOSUM80[i][j];
		break;
	}
	case 3:
	{
		int BLOSUM90[24][24] = {
		5,-2,-2,-3,-1,-1,-1,0,-2,-2,-2,-1,-2,-3,-1,1,0,-4,-3,-1,-2,-1,-1,-6,
		-2,6,-1,-3,-5,1,-1,-3,0,-4,-3,2,-2,-4,-3,-1,-2,-4,-3,-3,-2,0,-1,-6,
		-2,-1,7,1,-4,0,-1,-1,0,-4,-4,0,-3,-4,-3,0,0,-5,-3,-4,5,-1,-1,-6,
		-3,-3,1,7,-5,-1,1,-2,-2,-5,-5,-1,-4,-5,-3,-1,-2,-6,-4,-5,5,1,-1,-6,
		-1,-5,-4,-5,9,-4,-6,-4,-5,-2,-2,-4,-2,-3,-4,-2,-2,-4,-4,-2,-4,-5,-1,-6,
		-1,1,0,-1,-4,7,2,-3,1,-4,-3,1,0,-4,-2,-1,-1,-3,-3,-3,-1,5,-1,-6,
		-1,-1,-1,1,-6,2,6,-3,-1,-4,-4,0,-3,-5,-2,-1,-1,-5,-4,-3,1,5,-1,-6,
		0,-3,-1,-2,-4,-3,-3,6,-3,-5,-5,-2,-4,-5,-3,-1,-3,-4,-5,-5,-2,-3,-1,-6,
		-2,0,0,-2,-5,1,-1,-3,8,-4,-4,-1,-3,-2,-3,-2,-2,-3,1,-4,-1,0,-1,-6,
		-2,-4,-4,-5,-2,-4,-4,-5,-4,5,1,-4,1,-1,-4,-3,-1,-4,-2,3,-5,-4,-1,-6,
		-2,-3,-4,-5,-2,-3,-4,-5,-4,1,5,-3,2,0,-4,-3,-2,-3,-2,0,-5,-4,-1,-6,
		-1,2,0,-1,-4,1,0,-2,-1,-4,-3,6,-2,-4,-2,-1,-1,-5,-3,-3,-1,1,-1,-6,
		-2,-2,-3,-4,-2,0,-3,-4,-3,1,2,-2,7,-1,-3,-2,-1,-2,-2,0,-4,-2,-1,-6,
		-3,-4,-4,-5,-3,-4,-5,-5,-2,-1,0,-4,-1,7,-4,-3,-3,0,3,-2,-4,-4,-1,-6,
		-1,-3,-3,-3,-4,-2,-2,-3,-3,-4,-4,-2,-3,-4,8,-2,-2,-5,-4,-3,-3,-2,-1,-6,
		1,-1,0,-1,-2,-1,-1,-1,-2,-3,-3,-1,-2,-3,-2,5,1,-4,-3,-2,0,-1,-1,-6,
		0,-2,0,-2,-2,-1,-1,-3,-2,-1,-2,-1,-1,-3,-2,1,6,-4,-2,-1,-1,-1,-1,-6,
		-4,-4,-5,-6,-4,-3,-5,-4,-3,-4,-3,-5,-2,0,-5,-4,-4,11,2,-3,-6,-4,-1,-6,
		-3,-3,-3,-4,-4,-3,-4,-5,1,-2,-2,-3,-2,3,-4,-3,-2,2,8,-3,-4,-3,-1,-6,
		-1,-3,-4,-5,-2,-3,-3,-5,-4,3,0,-3,0,-2,-3,-2,-1,-3,-3,5,-4,-3,-1,-6,
		-2,-2,5,5,-4,-1,1,-2,-1,-5,-5,-1,-4,-4,-3,0,-1,-6,-4,-4,5,0,-1,-6,
		-1,0,-1,1,-5,5,5,-3,0,-4,-4,1,-2,-4,-2,-1,-1,-4,-3,-3,0,5,-1,-6,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-6,
		-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,1
		};
		i = FindMatrixIndex(a);
		j = FindMatrixIndex(b);
		return BLOSUM90[i][j];
		break;
	}
	case 4:
	{
		int PAM30[24][24] = {
		6,-7,-4,-3,-6,-4,-2,-2,-7,-5,-6,-7,-5,-8,-2,0,-1,-13,-8,-2,-3,-3,-1,-17,
		-7,8,-6,-10,-8,-2,-9,-9,-2,-5,-8,0,-4,-9,-4,-3,-6,-2,-10,-8,-7,-4,-1,-17,
		-4,-6,8,2,-11,-3,-2,-3,0,-5,-7,-1,-9,-9,-6,0,-2,-8,-4,-8,6,-3,-1,-17,
		-3,-10,2,8,-14,-2,2,-3,-4,-7,-12,-4,-11,-15,-8,-4,-5,-15,-11,-8,6,1,-1,-17,
		-6,-8,-11,-14,10,-14,-14,-9,-7,-6,-15,-14,-13,-13,-8,-3,-8,-15,-4,-6,-12,-14,-1,-17,
		-4,-2,-3,-2,-14,8,1,-7,1,-8,-5,-3,-4,-13,-3,-5,-5,-13,-12,-7,-3,6,-1,-17,
		-2,-9,-2,2,-14,1,8,-4,-5,-5,-9,-4,-7,-14,-5,-4,-6,-17,-8,-6,1,6,-1,-17,
		-2,-9,-3,-3,-9,-7,-4,6,-9,-11,-10,-7,-8,-9,-6,-2,-6,-15,-14,-5,-3,-5,-1,-17,
		-7,-2,0,-4,-7,1,-5,-9,9,-9,-6,-6,-10,-6,-4,-6,-7,-7,-3,-6,-1,-1,-1,-17,
		-5,-5,-5,-7,-6,-8,-5,-11,-9,8,-1,-6,-1,-2,-8,-7,-2,-14,-6,2,-6,-6,-1,-17,
		-6,-8,-7,-12,-15,-5,-9,-10,-6,-1,7,-8,1,-3,-7,-8,-7,-6,-7,-2,-9,-7,-1,-17,
		-7,0,-1,-4,-14,-3,-4,-7,-6,-6,-8,7,-2,-14,-6,-4,-3,-12,-9,-9,-2,-4,-1,-17,
		-5,-4,-9,-11,-13,-4,-7,-8,-10,-1,1,-2,11,-4,-8,-5,-4,-13,-11,-1,-10,-5,-1,-17,
		-8,-9,-9,-15,-13,-13,-14,-9,-6,-2,-3,-14,-4,9,-10,-6,-9,-4,2,-8,-10,-13,-1,-17,
		-2,-4,-6,-8,-8,-3,-5,-6,-4,-8,-7,-6,-8,-10,8,-2,-4,-14,-13,-6,-7,-4,-1,-17,
		0,-3,0,-4,-3,-5,-4,-2,-6,-7,-8,-4,-5,-6,-2,6,0,-5,-7,-6,-1,-5,-1,-17,
		-1,-6,-2,-5,-8,-5,-6,-6,-7,-2,-7,-3,-4,-9,-4,0,7,-13,-6,-3,-3,-6,-1,-17,
		-13,-2,-8,-15,-15,-13,-17,-15,-7,-14,-6,-12,-13,-4,-14,-5,-13,13,-5,-15,-10,-14,-1,-17,
		-8,-10,-4,-11,-4,-12,-8,-14,-3,-6,-7,-9,-11,2,-13,-7,-6,-5,10,-7,-6,-9,-1,-17,
		-2,-8,-8,-8,-6,-7,-6,-5,-6,2,-2,-9,-1,-8,-6,-6,-3,-15,-7,7,-8,-6,-1,-17,
		-3,-7,6,6,-12,-3,1,-3,-1,-6,-9,-2,-10,-10,-7,-1,-3,-10,-6,-8,6,0,-1,-17,
		-3,-4,-3,1,-14,6,6,-5,-1,-6,-7,-4,-5,-13,-4,-5,-6,-14,-9,-6,0,6,-1,-17,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-17,
		-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,1
		};
		i = FindMatrixIndex(a);
		j = FindMatrixIndex(b);
		return PAM30[i][j];
		break;
	}
	case 5:
	{
		int PAM70[24][24] = {
		5,-4,-2,-1,-4,-2,-1,0,-4,-2,-4,-4,-3,-6,0,1,1,-9,-5,-1,-1,-1,-1,-11,
		-4,8,-3,-6,-5,0,-5,-6,0,-3,-6,2,-2,-7,-2,-1,-4,0,-7,-5,-4,-2,-1,-11,
		-2,-3,6,3,-7,-1,0,-1,1,-3,-5,0,-5,-6,-3,1,0,-6,-3,-5,5,-1,-1,-11,
		-1,-6,3,6,-9,0,3,-1,-1,-5,-8,-2,-7,-10,-4,-1,-2,-10,-7,-5,5,2,-1,-11,
		-4,-5,-7,-9,9,-9,-9,-6,-5,-4,-10,-9,-9,-8,-5,-1,-5,-11,-2,-4,-8,-9,-1,-11,
		-2,0,-1,0,-9,7,2,-4,2,-5,-3,-1,-2,-9,-1,-3,-3,-8,-8,-4,-1,5,-1,-11,
		-1,-5,0,3,-9,2,6,-2,-2,-4,-6,-2,-4,-9,-3,-2,-3,-11,-6,-4,2,5,-1,-11,
		0,-6,-1,-1,-6,-4,-2,6,-6,-6,-7,-5,-6,-7,-3,0,-3,-10,-9,-3,-1,-3,-1,-11,
		-4,0,1,-1,-5,2,-2,-6,8,-6,-4,-3,-6,-4,-2,-3,-4,-5,-1,-4,0,1,-1,-11,
		-2,-3,-3,-5,-4,-5,-4,-6,-6,7,1,-4,1,0,-5,-4,-1,-9,-4,3,-4,-4,-1,-11,
		-4,-6,-5,-8,-10,-3,-6,-7,-4,1,6,-5,2,-1,-5,-6,-4,-4,-4,0,-6,-4,-1,-11,
		-4,2,0,-2,-9,-1,-2,-5,-3,-4,-5,6,0,-9,-4,-2,-1,-7,-7,-6,-1,-2,-1,-11,
		-3,-2,-5,-7,-9,-2,-4,-6,-6,1,2,0,10,-2,-5,-3,-2,-8,-7,0,-6,-3,-1,-11,
		-6,-7,-6,-10,-8,-9,-9,-7,-4,0,-1,-9,-2,8,-7,-4,-6,-2,4,-5,-7,-9,-1,-11,
		0,-2,-3,-4,-5,-1,-3,-3,-2,-5,-5,-4,-5,-7,7,0,-2,-9,-9,-3,-4,-2,-1,-11,
		1,-1,1,-1,-1,-3,-2,0,-3,-4,-6,-2,-3,-4,0,5,2,-3,-5,-3,0,-2,-1,-11,
		1,-4,0,-2,-5,-3,-3,-3,-4,-1,-4,-1,-2,-6,-2,2,6,-8,-4,-1,-1,-3,-1,-11,
		-9,0,-6,-10,-11,-8,-11,-10,-5,-9,-4,-7,-8,-2,-9,-3,-8,13,-3,-10,-7,-10,-1,-11,
		-5,-7,-3,-7,-2,-8,-6,-9,-1,-4,-4,-7,-7,4,-9,-5,-4,-3,9,-5,-4,-7,-1,-11,
		-1,-5,-5,-5,-4,-4,-4,-3,-4,3,0,-6,0,-5,-3,-3,-1,-10,-5,6,-5,-4,-1,-11,
		-1,-4,5,5,-8,-1,2,-1,0,-4,-6,-1,-6,-7,-4,0,-1,-7,-4,-5,5,1,-1,-11,
		-1,-2,-1,2,-9,5,5,-3,1,-4,-4,-2,-3,-9,-2,-2,-3,-10,-7,-4,1,5,-1,-11,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-11,
		-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,1
		};
		i = FindMatrixIndex(a);
		j = FindMatrixIndex(b);
		return PAM70[i][j];
		break;
	}
	case 6:
	{
		int PAM250[24][24] = {
		2,-2,0,0,-2,0,0,1,-1,-1,-2,-1,-1,-3,1,1,1,-6,-3,0,0,0,-1,-8,
		-2,6,0,-1,-4,1,-1,-3,2,-2,-3,3,0,-4,0,0,-1,2,-4,-2,-1,0,-1,-8,
		0,0,2,2,-4,1,1,0,2,-2,-3,1,-2,-3,0,1,0,-4,-2,-2,2,1,-1,-8,
		0,-1,2,4,-5,2,3,1,1,-2,-4,0,-3,-6,-1,0,0,-7,-4,-2,3,3,-1,-8,
		-2,-4,-4,-5,12,-5,-5,-3,-3,-2,-6,-5,-5,-4,-3,0,-2,-8,0,-2,-4,-5,-1,-8,
		0,1,1,2,-5,4,2,-1,3,-2,-2,1,-1,-5,0,-1,-1,-5,-4,-2,1,3,-1,-8,
		0,-1,1,3,-5,2,4,0,1,-2,-3,0,-2,-5,-1,0,0,-7,-4,-2,3,3,-1,-8,
		1,-3,0,1,-3,-1,0,5,-2,-3,-4,-2,-3,-5,0,1,0,-7,-5,-1,0,0,-1,-8,
		-1,2,2,1,-3,3,1,-2,6,-2,-2,0,-2,-2,0,-1,-1,-3,0,-2,1,2,-1,-8,
		-1,-2,-2,-2,-2,-2,-2,-3,-2,5,2,-2,2,1,-2,-1,0,-5,-1,4,-2,-2,-1,-8,
		-2,-3,-3,-4,-6,-2,-3,-4,-2,2,6,-3,4,2,-3,-3,-2,-2,-1,2,-3,-3,-1,-8,
		-1,3,1,0,-5,1,0,-2,0,-2,-3,5,0,-5,-1,0,0,-3,-4,-2,1,0,-1,-8,
		-1,0,-2,-3,-5,-1,-2,-3,-2,2,4,0,6,0,-2,-2,-1,-4,-2,2,-2,-2,-1,-8,
		-3,-4,-3,-6,-4,-5,-5,-5,-2,1,2,-5,0,9,-5,-3,-3,0,7,-1,-4,-5,-1,-8,
		1,0,0,-1,-3,0,-1,0,0,-2,-3,-1,-2,-5,6,1,0,-6,-5,-1,-1,0,-1,-8,
		1,0,1,0,0,-1,0,1,-1,-1,-3,0,-2,-3,1,2,1,-2,-3,-1,0,0,-1,-8,
		1,-1,0,0,-2,-1,0,0,-1,0,-2,0,-1,-3,0,1,3,-5,-3,0,0,-1,-1,-8,
		-6,2,-4,-7,-8,-5,-7,-7,-3,-5,-2,-3,-4,0,-6,-2,-5,17,0,-6,-5,-6,-1,-8,
		-3,-4,-2,-4,0,-4,-4,-5,0,-1,-1,-4,-2,7,-5,-3,-3,0,10,-2,-3,-4,-1,-8,
		0,-2,-2,-2,-2,-2,-2,-1,-2,4,2,-2,2,-1,-1,-1,0,-6,-2,4,-2,-2,-1,-8,
		0,-1,2,3,-4,1,3,0,1,-2,-3,1,-2,-4,-1,0,0,-5,-3,-2,3,2,-1,-8,
		0,0,1,3,-5,3,3,0,2,-2,-3,0,-2,-5,0,0,-1,-6,-4,-2,2,3,-1,-8,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-8,
		-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,1
		};
		i = FindMatrixIndex(a);
		j = FindMatrixIndex(b);
		return PAM250[i][j];
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

void ReverseFun(char* seqArr1, char* seqArr2, double simThr, int matrix, int exact)
{
	int row, col;
	int i = 0, j = 0;
	int n, a = 0;
	int k = 1, o = 1, l = 0;
	int count, similarityCount, identityCount;
	int diago, match;
	int getScore;

	double similarity, identity;
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
	similarityCount = 0;
	identityCount = 0;
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
				if (NonExactMatch(seqArr1[i - 1], revSeqArr[j - 1], matrix) >= 0)
				{
					while (k)
					{
						if ('F' == directs[i][j])
						{
							if (NonExactMatch(seqArr1[i - 1], revSeqArr[j - 1], matrix) >= 0 || ('X' == seqArr1[i - 1] && 'X' == revSeqArr[j - 1]))
							{
								similarityCount++;
							}
							if (seqArr1[i - 1] == revSeqArr[j - 1])
							{
								identityCount++;
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
							if (NonExactMatch(seqArr1[i + 1 - 1], revSeqArr[j + 1 - 1], matrix) >= 0)
							{
								similarity = (double)similarityCount / (double)l;
								identity = (double)identityCount / (double)l;
								if (similarity >= simThr && l >= 1)
								{
									fprintf(pfWrite, "%d\t%d\t%d\t%1f\t%d\t%1f\t%d\n", i+1, col-1-j, l, identity, l-identityCount, similarity, getScore);
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
									if (NonExactMatch(seqArr1[i + 1 - 1], revSeqArr[j + 1 - 1], matrix) >= 0)
									{
										similarity = (double)similarityCount / (double)l;
										identity = (double)identityCount / (double)l;
										if (similarity >= simThr && l >= 1)
										{
											fprintf(pfWrite, "%d\t%d\t%d\t%1f\t%d\t%1f\t%d\n", i+1, col-1-j, l, identity, l-identityCount, similarity, getScore);
										}
										break;
									}
								}
							}
							similarityCount = 0;
							identityCount = 0;
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

void DirectFun(char* seqArr1, char* seqArr2, double simThr, int matrix, int exact)
{
	int row, col;
	int i = 0, j = 0;
	int n, a = 0;
	int k = 1, o = 1, l = 0;
	int count, similarityCount, identityCount;
	int diago, match;
	int getScore;

	double similarity, identity;
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

	FILE* pfWrite = fopen("position.txt","w");
	count = 0;
	similarityCount = 0;
	identityCount = 0;
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
				if (NonExactMatch(seqArr1[i - 1], seqArr2[j - 1], matrix) >= 0)
				{
					while (k)
					{
						if ('F' == directs[i][j])
						{
							if (NonExactMatch(seqArr1[i - 1], seqArr2[j - 1], matrix) >= 0 || ('X' == seqArr1[i - 1] && 'X' == seqArr2[j - 1]))
							{
								similarityCount++;
							}
							if (seqArr1[i - 1] == seqArr2[j - 1])
							{
								identityCount++;
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
							if (NonExactMatch(seqArr1[i + 1 - 1], seqArr2[j + 1 - 1], matrix) >= 0)
							{
								similarity = (double)similarityCount / (double)l;
								identity = (double)identityCount / (double)l;
								if (similarity >= simThr && l >= 1)
								{
									fprintf(pfWrite, "%d\t%d\t%d\t%1f\t%d\t%1f\t%d\n", i+1, j+1, l, identity, l-identityCount, similarity, getScore);
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
									if (NonExactMatch(seqArr1[i + 1 - 1], seqArr2[j + 1 - 1], matrix) >= 0)
									{
										similarity = (double)similarityCount / (double)l;
										identity = (double)identityCount / (double)l;
										if (similarity >= simThr && l >= 1)
										{
											fprintf(pfWrite, "%d\t%d\t%d\t%1f\t%d\t%1f\t%d\n", i+1, j+1, l, identity, l-identityCount, similarity, getScore);
										}
										break;
									}
								}
							}
							similarityCount = 0;
							identityCount = 0;
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

	double simThr = atof(argv[1]);
	int mode = atoi(argv[2]);
	int matrix = atoi(argv[3]);
	int exact = atoi(argv[4]);

	if (0 == mode)
	{
		DirectFun(seqArr1, seqArr2, simThr, matrix, exact);
	}
	else if (1 == mode)
	{
		ReverseFun(seqArr1, seqArr2, simThr, matrix, exact);
	}


	free(seqArr1);
	free(seqArr2);

	return 0;

	// created by Huilong Chen!
}
