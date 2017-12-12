#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <sys/types.h>

double myclock();

__global__ 
void createMatrix(int** matrix, int rows, int columns, int *ksuCompact, int *sampleCompact, int *AB) {
	
	int i, j;	
		
	//printf("Creating matrix...\n");
	matrix[0][0] = 0;
	//printf("matrix[0][0] = %d\n", matrix[0][0]);
	for(i = 1; i < columns + 1;i++) {
		matrix[0][i] = ksuCompact[i - 1];
	}
	
	for(i = 1; i < rows + 1; i++) {
		matrix[i][0] = sampleCompact[i - 1];
		//printf("sampleCompact[%d] = %d, matrix[%d][0] = %d\n", i-1, sampleCompact[i-1], i, matrix[i][0]);	
	}

	for(j = 1; j < columns + 1; j++){
		for(i = 1;i < rows + 1;i++) {
			matrix[i][j] = AB[(i - 1) + (j - 1)*(rows)];
			//printf("matrix[%d][%d] = %d, AB[%d] = %d\n", i, j, matrix[i][j], ((i - 1) + (j - 1)*rows), AB[(i - 1) + (j - 1)*rows]);
		} 
	}
	printf("finished createMatrix.\n");
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
	
//	double tstart;
//	double ttotal;
//	struct rusage r_usage;	

	int i, j, err;

	FILE *fd;
	
	int maxlines = atoi(argv[1]);
	int nlines;
	
	//nColumns and nRows refers to the number of columns and rows of DATA not total in matrix.
	int nColumns;
	int nRows;
	
	int** matrixAB;

	int *ksuID;
	int *sampleID;
	int *genotypeAB;

	int *ksuIDCompact;
	int *sampleIDCompact;

//	tstart = myclock();
//	tstart = myclock();

	printf("Allocating memory...\n");
	fflush(stdout);
	
	//allocate memory for each line
		
	ksuID = (int*) malloc(maxlines * sizeof(int));
	sampleID = (int*) malloc(maxlines * sizeof(int));
	genotypeAB = (int*) malloc(maxlines * sizeof(int));
	
	ksuIDCompact = (int*) malloc(maxlines * sizeof(int));
	sampleIDCompact = (int*) malloc(maxlines * sizeof(int));
	
	/*	
	cudaMallocManaged(&ksuID, maxlines * sizeof(int));
	cudaMallocManaged(&sampleID, maxlines * sizeof(int));
	cudaMallocManaged(&genotypeAB, maxlines * sizeof(int));

	cudaMallocManaged(&ksuIDCompact, maxlines * sizeof(int));
	cudaMallocManaged(&sampleIDCompact, maxlines * sizeof(int));
	*/

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
	
	//matrix = (int**) malloc((rows+1) * sizeof(int*));
	cudaMallocManaged(&matrixAB, (nRows+1) * sizeof(int*));
	
	for(i = 0; i < nRows+1; i++) 
		//matrix[i] = (int*) malloc((columns+1) * sizeof(int));
		cudaMallocManaged(&matrixAB[i], (nColumns+1) * sizeof(int));

	cuda_error_t	
	createMatrix<<<1,1>>>(matrixAB, nRows, nColumns, ksuIDCompact, sampleIDCompact, genotypeAB);
	cudaDeviceSynchronize();

	printf("Matrix created.\n");

	//Write matrix to file
	fd = fopen("./SNPmatrix.txt", "w");
	if (fd != NULL) {
		printf("SNPmatrix.txt created/opened\n");
		//create first row
		for(i = 0; i < nRows + 1; i++) {
			for(j = 0; j < nColumns; j++) {
				printf("i = %d, j = %d\n, matrix[i][j] = %d\n", i, j, matrixAB[i][j]);
				fprintf(fd, "%d\t", matrixAB[i][j]);
			}
			fprintf(fd, "%d\n", matrixAB[i][nColumns]);
		}		
		fclose(fd);
	}
	printf("Closed SNPmatrix.txt.\n");
//	ttotal = myclock() - tstart;
	//getrusage(RUSAGE_SELF, &r_usage);
	
//	printf("SNPmatrix.txt created.\nExec. Time = %f, RAM Usage = %ld", ttotal, r_usage.ru_maxrss);
		

	free(ksuID);
	free(sampleID);
	free(genotypeAB);
	free(sampleIDCompact);
	free(ksuIDCompact);
	cudaFree(matrixAB);
	
	return 1;
}

double myclock() {
	static time_t t_start = 0;

	struct timespec ts;
	clock_gettime(CLOCK_REALTIME, &ts);
	if( t_start == 0 ) t_start = ts.tv_sec;

	return (double) (ts.tv_sec - t_start) + ts.tv_nsec * 1.0e-9;
}
