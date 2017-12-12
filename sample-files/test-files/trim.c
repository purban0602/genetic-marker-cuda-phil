#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <sys/resource.h>

double myclock();

int main(int argc, char **argv) {
	
	int myID, nProcs;
	int nCores, nNodes;

	int nwords, maxwords = 50000;
	int nlines, maxlines = 1000000;
	
	int i, k, n, err, done, *count;
	
	//Each process spins, requesting the next baatch of keywords from rank 0.
	//int myBlockStart, myBlockEnd, blockSize;
	int myBlockStart;
	int blockSize = 100;

	double nchars = 0;
	double tstart, ttotal;
	FILE *fd;

	char **word, **line;
	
	MPI_Status status;

	maxlines = atoi(argv[1]);

	struct rusage r_usage;

	err = MPI_Init(&argc, &argv);

	err = MPI_Comm_rank(MPI_COMM_WORLD, &myID);
	err = MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

	if (myID == 0) {

		tstart = myclock();
		tstart = myclock();
		//Allocate memory for each line in file
		line = (char **) malloc(maxlines * sizeof(char*));
	
		for(i = 0; i < maxlines; i++) {
			line[i] = malloc(2001);
		}
	
		//open genotypes file
		fd = fopen("/home/phillip/Documents/Genotypes6.txt", "r");
		if(fd != NULL) {
			printf("Opened Genotypes6.txt\n");
			fflush(stdout);
			//Store lines in line array
			nlines = -1;
			do {
				err = fscanf(fd, "%[^\n]\n", line[++nlines]);
				if(line[nlines] != NULL) nchars += (double) strlen(line[nlines]);
			} while(err != EOF && nlines < maxlines);
				
			printf("Closing Genotypes6.txt...\n");
			fflush(stdout);
			fclose(fd);
		}
		printf("Now writing rawdata.txt. line[0] = \"%s\", nwords = %d\n", line[0], nlines);
		fflush(stdout);	
		fd = fopen("/home/phillip/Documents/genetic-marker-cuda/sample-files/test-files/rawdata.txt", "w");
		for(k = 0; k < maxlines && k < nlines; k++) {
			fprintf(fd, "%s\n", line[k]);
		}
		fclose(fd);
		printf("Wrote rawdata.txt.\n");
		fflush(stdout);
	
		ttotal = myclock() - tstart;

		getrusage(RUSAGE_SELF, &r_usage);
		printf("%f, %ld\n", ttotal, r_usage.ru_maxrss);
	}
	
	MPI_Finalize();
}

double myclock() {
	static time_t t_start = 0;

	struct timespec ts;
	clock_gettime(CLOCK_REALTIME, &ts);
	if( t_start == 0) t_start = ts.tv_sec;

	return (double) (ts.tv_sec - t_start) + ts.tv_nsec * 1.0e-9;
}                                                                                                                    
