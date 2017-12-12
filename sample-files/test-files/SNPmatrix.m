

% Input numeric file (header row trimmed) wtih KSUid, SampleID, and
% genotypeAB in columns 1-3 respectively and title RawData.  Output is SNP data matrix with
% animals in rows and SNPs in columns.
fileID = fopen('rawdata.txt');
RawData = textscan(fileID,'%sq');

NumAnimals=size(unique(RawData(:,2)),1);
NumSNPs=size(unique(RawData(:,1)),1);
DataMatrix=zeros(NumAnimals, NumSNPs);

SNPs=unique(RawData(:,1));
Animals=unique(RawData(:,2));
DataMatrix=[Animals DataMatrix];
Temp=[0; SNPs]';
DataMatrix=[Temp; DataMatrix];
clear Temp

%Typically I would sort and dump genotypes in as rows so that it doens't take ages to complete, 
%but not sure if you
%want to do it that way (espeically if you want to parallelize), so I'm
%programming it out the long way in case there's a better solution

for a=1:NumAnimals;
    AnimalID=Animals(a,1);
    Row=find(DataMatrix(:,1)==AnimalID);
    for b=1:NumSNPs;
        SNPid=SNPs(b,1);
        Geno=RawData(find(RawData(:,1)==SNPid & RawData(:,2)==AnimalID),3);
        Col=find(DataMatrix(1,:)==SNPid);
        DataMatrix(Row,Col)=Geno;
    end
end
clear Row
clear Col
clear Geno
clear AnimalID
clear SNPid
        
        
        
        



