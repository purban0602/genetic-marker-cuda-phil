#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <sys/types.h>

double myclock();


void main(int argc, char **argv) {
	
	int numAnimals, numSNPs = 0;
	int i, j, k, err;

	FILE *fd;
	
	int maxlines = atoi(argv[1]);
	int nlines;
	int nColumns;
	int nRows;
	int matrixN;

	int *ksuID;
	int *sampleID;
	int *genotypeAB;

	char **line;
	char tempLine[50];

	struct rusage r_usage;
	printf("Allocating memory...\n");
	fflush(stdout);
	//allocate memory for each line
	
	ksuID = (int*) malloc(maxlines * sizeof(int*));
	sampleID = (int*) malloc(maxlines * sizeof(int*));
	genotypeAB = (int*) malloc(maxlines * sizeof(int*));

	//assume input file is already 3 columns needed for data matrix
	/*

	//Create "rawdata" file from "FinalReport" file.
	//open raw data file (KSUid, SampleID, genotypeAB)
	printf("Opening FinalReport_Truncated.txt...\n");
	fflush(stdout);

	fd = fopen("./FinalReport_Truncated.txt","r");
	
	printf("Opened FinalReport_Truncated.txt.\nRemoving header...\n");
	fflush(stdout);
	
	if (fd != NULL) {
		while (strcmp("Y",tempLine) != 0) {
			err = fscanf(fd, "%s", tempLine);	
			printf("Current string is %s.\n", tempLine);

		}
		printf("Header removed.\nConverting data...\n");

		nlines = -1;
		do {
			err = fscanf(fd, "%s\t%d\t%c\t%c\t%c\t%c\t%c\t%c\t%f\t%f\t%f\n", tempLine);
			printf("Current string is %s.\n", tempLine);
			nlines++;
		} while(err != EOF && nlines < maxlines);
		fclose(fd);
	}
	else {
		printf("Error: File Not Found.\n");
	}
	*/
	
	//Only need the unique animals and SNPs
	fd = fopen("./rawdata.txt", "r");
	printf("Opened rawdata.txt.\n");	
	if (fd != NULL) {
		nlines = 0;
		do {
			err = fscanf(fd, "%d", &ksuID[nlines]);		
			err = fscanf(fd, "%d", &sampleID[nlines]);		
			err = fscanf(fd, "%d", &genotypeAB[nlines]);

			printf("%d %d %d\n", ksuID[nlines], sampleID[nlines], genotypeAB[nlines]);
			nlines++;
		} while(err != EOF && nlines < maxlines);
		fclose(fd);
	}
	printf("File read successfully.\nWriting matrix...\n");
	fflush(stdout);
	/* output matrix: ksu ids in columns, sample ids in rows, genotypeAB everywhere else
	 * matrix looks like this...
	 * 	0	1	2	3	...
	 * 	737	2	1	2	...
	 *	926	2	1	3	...
	 *	948	3	1	1	...
	 *	...	...	...	...
	 */	

	fd = fopen("./testmatrix.txt", "w");

	//create first row
	fprintf(fd, "0\t");
	fprintf(fd, "%d\t",ksuID[0]);
	nColumns = 1;
	
	for (i = 1; i < nlines-1; i++) {
		if (ksuID[i] != ksuID[i-1]) {
			fprintf(fd, "%d\t",ksuID[i]);
			nColumns++;
		}
	}

	//avoids an extra \t at the end of the row
	if (ksuID[nlines-1] != ksuID[nlines-2]) {
	       fprintf(fd, "%d",ksuID[i]);  
	       nColumns++;
	}

	fprintf(fd, "\n");	
	printf("nColumns = %d\n", nColumns);
	
	if(nlines%nColumns) nRows = (nlines / nColumns) + 1;
	else nRows = nlines/nColumns;
	printf("nRows = %d\n", nRows);
	printf("nlines = %d\n", nlines);

	//create all other rows
	matrixN = 0;
	for (i = 0; i < nRows; i++) {
		//write first column (sampleID)
		fprintf(fd, "%d\t", sampleID[i]);

		//write all other columns
		for (j = 0; j < nColumns-1; j++) {
			fprintf(fd, "%d\t", genotypeAB[i + (j*nRows)]);
			matrixN++;
		}
		fprintf(fd, "%d\n", genotypeAB[i + ((nColumns - 1)*nRows)]); 
		
	}
	
}
