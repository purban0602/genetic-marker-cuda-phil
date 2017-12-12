Caleb, I've tried to provide some notes on workflow and what needs to be done 
in what order-please let me know if something doens't make sense!  I realize you may not be
doing all of these steps, but this is what I have envisioned in my head and what I hope to 
achieve at some point-whatever you guys have agreed to focus on for your project is perfectly fine.  

1.  In MapFile, recode chromsome values as follows: X recode as 30, Y recode as 31 and 0 recode as 32
	Reorder entries by sorting first by chromosome and then by position. 
 	Create KSUid in column 1 using genomic position to order
	Create column at the end (CumGenomicPos)-first entry equal to position column, next rows equal to the CumGenomicPos of entry above plus Position in same row
	Write out FinalMapFile.txt
2.  In Final Report, replace SNP Name column for each row with corresponding KSUid from MapFile
	Create new column for genotypeAB
	Recode genotypes as follows:  	if Allele1AB=A and Allele2AB=B, genotypeAB=3
					if Allele1AB=A and Allele2AB=A, genotypeAB=1
					if Allele1AB=B and Allele2AB=B, genotypeAB=2
					if Allele1AB=B and Allele2AB=A, genotypeAB=3
					If Allele1AB and/or Allele2AB are missing, genotypeAB=99 (both should be missing if one is missing due to way genotypes are called)
	Write out this file FinalReport_Numeric.txt
3.  On Final Report, delete header rows (9 rows, 10th is headers for the actual data)
	Delete columns SNP Name (if replaced by KSUid, leave in file), Allele1 Forward, Allele2 Top, Allele1 Top, Allele2 Top
	Delete columns Allele1AB, Allele2AB, GC Score, X, and Y
	Remaining columns should be KSUid, SampleID, genotypeAB

4.  Create data matrix with rows as animals and SNPs as columns. (this is my Matlab code-SNPmatrix.m)
	Place the genotypeAB from the FinalReport for each animal and SNP in the appropriate columns
	Write out data matrix as SNPmatrix.txt (DataMatrix_Truncated.txt in the sample files for comparison)
