#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <sys/types.h>

double myclock();

__global__ 
void createMatrix(int** matrix, int rows, int columns, int *ksuCompact, int *sampleCompact) {
	int i, j;
	//fprintf(fd, "0\t");
	
	//matrix = (int**) malloc((rows+1) * sizeof(int*));
	cudaMallocManaged(&matrix, (rows + 1) * sizeof(int*));
	
	
	for(i = 0; i < rows; i++) 
		//matrix[i] = (int*) malloc((columns+1) * sizeof(int));
		cudaMallocManaged(&matrix[i], (columns+1) * sizeof(int*));
	
	matrix[0][0] = 0;
	printf("matrix[0][0] = %d\n", matrix[0][0]);
}

int getColumns(int *ksu, int lines, int *ksuCompact) {

	int i, j;
	int exists;
	int compactIndex = 1;
	
	//create compact array (unique ksuID's)
	//place first element
	ksuCompact[0] = ksu[0];
//	printf("ksuCompact[0] = %d\n", ksu[0]);
	for(i = 1; i < lines; i++) {
		exists = 0;
		for(j = 0; j < i; j++) {
//			printf("ksuCompact[j] = %d, ksu[i] = %d\n", ksu[j], ksu[i]);
			if(ksu[i] == ksuCompact[j]) {
				exists = 1;
			}
		}
		if (!exists) {
			ksuCompact[compactIndex] = ksu[i];
			compactIndex++;
		}
	}

	return compactIndex;
}

int getRows(int lines, int columns, int *sample, int *sampleCompact) {

	int i, j;
	int exists;
	int compactIndex = 1;

	//create compact array (unique sampleID's)
	//place first element
	sampleCompact[0] = sample[0];
//	printf("sampleCompact[0] = %d\n", sample[0]);
	for(i = 1; i < lines; i++) {
		exists = 0;
		for(j = 0; j < i; j++) {
//			printf("sampleCompact[j] = %d, sample[i] = %d\n", sampleCompact[j], sample[i]);
			//if sample[i] currently exists in sampleCompact
			if(sample[i] == sampleCompact[j]) {
				exists = 1;
			}
		}
		if (!exists) {
			sampleCompact[compactIndex] = sample[i];
			compactIndex++;
		}
	}

	return compactIndex;

}

int main(int argc, char **argv) {
	
	int numAnimals, numSNPs = 0;
	int i, j, k, err;

	FILE *fd;
	
	int maxlines = atoi(argv[1]);
	int nlines;

	//nColumns and nRows refers to the number of columns and rows of DATA not total in matrix.
	int nColumns;
	int nRows;
	int matrixN;

	int** matrixAB;

	int *ksuID;
	int *sampleID;
	int *genotypeAB;

	int *ksuIDCompact;
	int *sampleIDCompact;

	char **line;
	char tempLine[50];

	struct rusage r_usage;
	printf("Allocating memory...\n");
	fflush(stdout);
	//allocate memory for each line
	/*	
	ksuID = (int*) malloc(maxlines * sizeof(int*));
	sampleID = (int*) malloc(maxlines * sizeof(int*));
	genotypeAB = (int*) malloc(maxlines * sizeof(int*));
	
	ksuIDCompact = (int*) malloc(maxlines * sizeof(int*));
	sampleIDCompact = (int*) malloc(maxlines * sizeof(int*));
	*/

	cudaMallocManaged(&ksuID, maxlines * sizeof(int*));
	cudaMallocManaged(&sampleID, maxlines * sizeof(int*));
	cudaMallocManaged(&genotypeAB, maxlines * sizeof(int*));

	cudaMallocManaged(&ksuIDCompact, maxlines * sizeof(int*));
	cudaMallocManaged(&sampleIDCompact, maxlines * sizeof(int*));

	//assume input file is already 3 columns needed for data matrix
	
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

	//get number of columns
	nColumns = getColumns(ksuID, nlines, ksuIDCompact);
	printf("nColumns = %d\n", nColumns);
	//get number of rows
	nRows = getRows(nlines, nColumns, sampleID, sampleIDCompact);
	printf("nRows = %d\n", nRows);
	
	createMatrix<<<1, 1>>>(matrixAB, nRows, nColumns, ksuIDCompact, sampleIDCompact);
	//Write matrix to file
	//fd = fopen("./SNPmatrix.txt", "w");

	//create first row
	/*
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
	*/

	//printf("nRows = %d\n", nRows);
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
	free(ksuID);
	free(sampleID);
	free(genotypeAB);
	free(sampleIDCompact);
	free(ksuIDCompact);
	free(matrixAB);
	
	return 1;
}
