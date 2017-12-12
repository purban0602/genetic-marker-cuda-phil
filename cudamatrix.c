#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <sys/types.h>

double myclock();
int err;

/* output matrix: ksu ids in columns, sample ids in rows, genotypeAB everywhere else
 * matrix looks like this...
 * 	0	1	2	3	...
 * 	737	2	1	2	...
 *	926	2	1	3	...
 *	948	3	1	1	...
 *	...	...	...	...
 */	

void initMatrix(int **matrix, int rows, int columns) {
	int i;
	int j;

	for(i = 0; i < rows; i++)
		for(j = 0; j < columns; j++)
			matrix[i][j] = 99;

}

void fillMatrix(int **matrix, int lines, int *ksu, int *sample, int *AB, int rows, int columns) {
	
	//each index has a ksuID, sampleID, and genotypeAB
	//nlines is the amount of data
	//
	//for each index, check if ksuid exists already,
	//then, check if sample id exists already,
	//finally, put AB in the correct column
	int i, j, k;
	
	j = 1;
	//fills first row, assumes already sorted
	for(i = 1; i < lines; i++) {
		printf("trying i = %d, j = %d\n", i, j);
		if(ksu[i] != matrix[0][j-1]) {
			matrix[0][j] = ksu[i];
			printf("Column %d = %d\n", j, matrix[0][j]);
			j++;
		}
	}

	printf("First row and first column are initialized.\n");
	/*
	//for each element...
	for(i = 0; i < lines; i++) {
		//check each column
		for(j = 0; j < columns; j++) {
			if(ksu[i] == matrix[0][j]) {
				//check each row
				for (k = 1; k < rows; k++) {
					if(sample[i] == matrix[k][0]) {
				//		matrix[k][j] = AB[i]
					}
				}
			}
		}
	}*/	
}

//currently not using!
int* createMatrix(int lines, int *ksu, int *sample, int *AB, int *pr, int *pc) {

	//initialize matrix
	int i, j;
	int nRows;
	int nColumns = 2;
	int *matrix;

	for (i = 1; i < lines; i++) {
		if(ksu[i] != ksu[i-1]) {
			nColumns++;
		}
	}

	nRows = (lines/(nColumns-1) + 1);
	//if there are empty elements
	if(lines%(nColumns-1) != 0) nRows++;

	matrix = (int*) malloc(nRows * nColumns * sizeof(int));	
	
	*pr = nRows;
	*pc = nColumns;
	printf("nRows = %d\nnColumns = %d\n", *pr, *pc);

	return matrix;
}

//reads input file, returns the number of lines in file
int readData(FILE *f, int lines, int max, int *ksu, int *sample, int *AB) {

	int i, j;

	//Only need the unique animals and SNPs
	f = fopen("./rawdata.txt", "r");
	printf("Opened rawdata.txt.\n");	
	if (f != NULL) {
		lines = 0;
		do {
			err = fscanf(f, "%d", &ksu[lines]);		
			err = fscanf(f, "%d", &sample[lines]);		
			err = fscanf(f, "%d", &AB[lines]);

			printf("%d %d %d\n", ksu[lines], sample[lines], AB[lines]);
			lines++;
		} while(err != EOF && lines < max);
		fclose(f);
	}
	return lines;
}


void main(int argc, char **argv) {
	
	int numAnimals, numSNPs = 0;
	int i, j;

	FILE *fd;
	
	int maxlines = atoi(argv[1]);
	int nlines;
	int nColumns;
	int nRows;
	int matrixN;

	int *prow = &nRows;
	int *pcolumn = &nColumns;

	int *ksuID;
	int *sampleID;
	int *genotypeAB;
	
	//matrix is 1D, must address as 2D
	int **matrixAB;

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
	
	//try to read data file
	nlines = readData(fd, nlines, maxlines, ksuID, sampleID, genotypeAB);

	printf("File read successfully.\nnlines = %d\nWriting matrix...\n", nlines);

	//initialize and create matrix...
	
	//matrixAB = createMatrix(nlines, ksuID, sampleID, genotypeAB, prow, pcolumn);
	
	//initialize matrix
	nColumns = 2;

	for (i = 1; i < nlines; i++) {
		if(ksuID[i] != ksuID[i-1]) {
			nColumns++;
		}
	}
	if(nlines%(nColumns-1)) nRows = nlines/((nColumns-1) + 2);
	else nRows = nlines/(nColumns-1) + 1;

	matrixAB = (int**) malloc(nRows * sizeof(int*));
	for(i = 0; i < nRows; i++) matrixAB[i] = (int*) malloc(nColumns * sizeof(int));

	printf("nRows = %d\nnColumns = %d\n", *prow, *pcolumn);
	
	initMatrix(matrixAB, nRows, nColumns);

	matrixAB[0][0] = 0;
	printf("matrixAB[0][0] = %d\n", matrixAB[0][0]);
	printf("matrixAB[1][1] = %d\n", matrixAB[1][1]);
	
	fillMatrix(matrixAB, nlines, ksuID, sampleID, genotypeAB, nRows, nColumns);
	
	/* output matrix: ksu ids in columns, sample ids in rows, genotypeAB everywhere else
	 * matrix looks like this...
	 * 	0	1	2	3	...
	 * 	737	2	1	2	...
	 *	926	2	1	3	...
	 *	948	3	1	1	...
	 *	...	...	...	...
	 */	
	
	
	/*
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
	*/
	
}
