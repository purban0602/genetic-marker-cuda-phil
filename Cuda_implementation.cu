#include <cuda_runtime.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/resource.h>


typedef struct{
	char** name;
	char* chrom_c;
	//int* chrom;
	long* pos;
	//long* c_pos;
	char** rest;
}SNP;

typedef struct{
	char* snp_name;
	int* a_id; //length is the number of animals
	char* ab1; 
	char* ab2;
	int* ab;
}Sample;

int NSNPS;
int NSAMPLES;

__device__ void sort_by_bit(SNP* snps, Sample* samples, int bit);

__device__ long scan(long* x);



void read_files(char* map_path, char* snp_path, char** data_string, char** snps_data){
	
	FILE *fd;
	int err;
	int num_lines = -1;
	char** header_array;
	int i;
	
/***********************Allocate string for header info**********/
	printf("Allocating string for header array...\n");

	header_array = (char**) malloc( 10 * sizeof(char*));
	
	for(i = 0; i < 10; i++){
		header_array[i] = (char*)malloc(100); 	
	}
/*****************************************************************/
	
	fd = fopen(snp_path, "r");
	
	
/*******Getting number of SNP and Sample from header****/
	printf("Getting number of SNPs and Samples from header...\n");

	do {
		err = fscanf(fd, "%[^\n]\n", header_array[++num_lines]);
	} while(err != EOF && num_lines < 10);
	
	
	err = sscanf(header_array[5], "Total SNP	%d", &NSNPS);
	err = sscanf(header_array[7], "Total Sample	%d", &NSAMPLES);
/***********************************************************/
	

	
/*************Getting Final Report Data***********************************/
	printf("Getting final report data...\n");

	//char** data_string;
	
	data_string = (char**) malloc(NSNPS * NSAMPLES * sizeof(char*));
	for(i = 0; i < NSNPS*NSAMPLES; i++){
		data_string[i] = (char*)malloc(100); 	
	}
	
	num_lines  = -1;
	do {
		err = fscanf(fd, "%[^\n]\n", data_string[++num_lines]);
	} while(err != EOF && num_lines < NSNPS*NSAMPLES);
	
	fclose(fd);
/**************************************************************************/

	
/************************Getting MapFile Data******************************/
	printf("Getting mapfile data...\n");

	//char** snps_data;
	char* junk = (char*) malloc(50 * sizeof(char));
	
	snps_data = (char**) malloc(NSNPS * sizeof(char*));
	for(i = 0; i < NSNPS; i++){
		snps_data[i] = (char*)malloc(100); 	
	}
	
	fd = fopen(map_path, "r");
	
	int num_lines2 = -1;
	err = fscanf(fd, "%[^\n]\n", junk);
	do {
		err = fscanf(fd, "%[^\n]\n", snps_data[++num_lines2]);
	} while(err != EOF && num_lines2 < NSNPS);
	
	free(junk);
	
	fclose(fd);
/**************************************************************************/
	

}

/*************functions for the radix sort**********************************/

__device__ void radixsort(SNP* snps, Sample* samples){
	
	for(int i = 0; i < 64; i++){
		sort_by_bit(snps, samples, i);
		__syncthreads();
	}
	
}

__device__ void sort_by_bit(SNP* snps, Sample* samples, int bit){
	
		int i = threadIdx.x;
		int size = blockDim.x;
		int index;
		
		/***temperary variables for the snps*****/
		long t_pos = snps->pos[i];
		char* t_name = snps->name[i];
		char t_chrom_c = snps->chrom_c[i];
		//char* t_rest = snps->rest[i];
		Sample t_sample = samples[i];
		
		int p_i = (t_pos >> bit) & 1;
		
		snps->pos[i] = p_i;
		
		__syncthreads();
		
		int ones_before = scan(snps->pos);
		int ones_total = snps->pos[size -1];
		int zeros_total = size - ones_total;
		
		__syncthreads();
		
		if(p_i)
			index = ones_before - 1 + zeros_total;
		else
			index = i - ones_before;
		
		snps->pos[index] = t_pos;
		snps->name[index] = t_name;
		snps->chrom_c[index] = t_chrom_c;
		//snps->rest[index] = t_rest;
		samples[index] = t_sample;
}

/**************************************************************************/

__device__ long scan(long* x){
	
	int i = threadIdx.x;
	int n = blockDim.x;
	int offset;
	
	for ( offset = 1; offset < n; offset *= 2){
		long temp;
		if (i >= offset)
			temp = x[i-offset];
		
		__syncthreads();
		
		if(i >= offset)
			x[i] = temp + x[i];
		
		__syncthreads();
	}
	
	return x[i];
}


void parse(SNP* snps, Sample* animals, char** data_string, char** snp_data){
	
	int i, j, err;
	
	snps->name = (char**) malloc(NSNPS * sizeof(char*));
	snps->chrom_c = (char*) malloc(NSNPS * sizeof(char));
	snps->pos = (long*) malloc(NSNPS * sizeof(long));
	
	for(i = 0; i < NSNPS; i++)
		snps->name[i] = (char*) malloc(50 * sizeof(char));
	
	animals = (Sample*) malloc(NSNPS * sizeof(Sample));
	
	for(i = 0; i < NSNPS; i++){
		animals[i].snp_name = (char*) malloc(50 * sizeof(char));
		animals[i].a_id = (int*) malloc(NSAMPLES * sizeof(int));
		animals[i].ab1 = (char*) malloc(NSAMPLES * sizeof(char));
		animals[i].ab2 = (char*) malloc(NSAMPLES * sizeof(char));
		animals[i].ab = (int*) malloc(NSAMPLES * sizeof(char));
	}
	
	for (i = 0; i < NSNPS; i++){
		err = sscanf(snp_data[i], "%*d	%s	%c	%ld	%*s", 
					  snps->name[i], snps->chrom_c[i], snps->pos[i], snps->rest[i]);
	}
	
	for(i = 0; i < NSNPS; i++){
		for(j = 0; j < NSAMPLES; j++)
			err = sscanf(data_string[i], "%s/t%d/t%*c/t%*c/t%*c/t%*c/t%c/t%c/t%*s", 
							animals[i].snp_name, animals[i].a_id[j], animals[i].ab1[j], animals[i].ab2[j]);
	}
}

__global__ void sort(SNP* snps, Sample* samples, int nsamples){
	
	int id = threadIdx.x;
	radixsort(snps, samples);
	
	for(int i = 0; i < nsamples; i++){
		if (samples[id].ab1[i] == 'A' && samples[id].ab2[i] == 'A'){
			samples[id].ab[i] = 1;
		}else if(samples[id].ab1[i] == 'B' && samples[id].ab2[i] == 'B'){
			samples[id].ab[i] = 2;
		}else{
			samples[id].ab[i] = 3;
		}
	}
}
int main(int argc, char** argv){
	printf("Begin.\n");
	
	SNP h_snps;
/*	
	typedef struct{
		char** name;
		char* chrom_c;
		//int* chrom;
		long* pos;
		//long* c_pos;
		//char** rest;
	}SNP;
*/
	Sample* h_samples;
/*
	typedef struct{
		char* snp_name;
		int* a_id; //length is the number of animals
		char* ab1; 
		char* ab2;
		int* ab;
	}Sample;
*/
	//char map_path[], snp_path[];
	char** data_string, **snps_data;
	char** d_name;
	char* d_chrom_c;
	long* d_pos;

	printf("Reading files...\n");
	
	//map_path = argv[1];
	char map_path[] = "./sample-files/test-files/SNP_Map_Truncated.txt";	
	//snp_path = argv[2];
	char snp_path[] = "./sample-files/test-files/FinalReport_Truncated.txt";
	
	read_files(map_path, snp_path, data_string, snps_data);

	printf("Files read.\nParsing...\n");

	parse(&h_snps, h_samples, data_string, snps_data);

	printf("Data parsed.\n");
	
	free(data_string);
	free(snps_data);
	
	printf("Allocating CUDA memory...\n");

	cudaMalloc((void**)&(d_pos), sizeof(long)*NSNPS);
	cudaMalloc((void**)&(d_chrom_c), sizeof(char)*NSNPS);
	cudaMalloc((void**)d_name, sizeof(char*)*NSNPS);
	
	cudaMemcpy(d_pos, (h_snps.pos), sizeof(long)*NSNPS, cudaMemcpyHostToDevice);
	cudaMemcpy(d_chrom_c, (h_snps.chrom_c), sizeof(char)*NSNPS, cudaMemcpyHostToDevice);
	cudaMemcpy(d_chrom_c, (h_snps.chrom_c), sizeof(char)*NSNPS, cudaMemcpyHostToDevice);
	
}
