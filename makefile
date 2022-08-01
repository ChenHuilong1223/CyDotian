CC = gcc

all: bpRepeatScan aaRepeatScan

bpRepeatScan: bpRepeatScan.c
	$(CC) -g -std=c99 bpRepeatScan.c -o bpRepeatScan

aaRepeatScan: aaRepeatScan.c
	$(CC) -g -std=c99 aaRepeatScan.c -o aaRepeatScan
