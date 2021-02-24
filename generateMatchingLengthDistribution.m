function matchingDataset = generateMatchingLengthDistribution(originalDataset,refDataset,sampleSize)

%% Martin White December 2020

% REF: XXXX

%Input: use the cell arrays created by generateComponentHomologs.m
%%


%Step 1: calcualte PDF (MLE gamma) of component homolog length 
%distributions in the reference dataset
refLengthGamFit = gamfit(refDataset{1,1}{1,1}(:,1));

%Step 2: For each component homolog in the dataset of interest, calculate 
%the probability of observing that length in the reference dataset

totalCHs = length(originalDataset{1,1}{1,1}(:,1));
probBasedOnRef(1:totalCHs,1)=nan;
for i = 1:totalCHs
    probBasedOnRef(i,1) = gampdf(originalDataset{1,1}{1,1}(i,1),refLengthGamFit(1),refLengthGamFit(2));
end


%Step 3: sample the original dataset of interest in a manner to best match
%the component homolog length distribution in the reference dataset
sampleIndex = randsample(1:length(originalDataset{1,1}{1,1}(:,1)),sampleSize,true,probBasedOnRef(:,1));

matchingDataset{1,1}=[];
matchingDataset{2,1}=[];
matchingDataset{3,1}=[];
matchingDataset{4,1}=[];
matchingDataset{5,1}=[];
matchingDataset{6,1}=[];

for i = 1:sampleSize
    
    matchingDataset{1,1}(i,:) = originalDataset{1,1}{1,1}(sampleIndex(i),:);
    matchingDataset{2,1}(i,:) = originalDataset{1,1}{2,1}(sampleIndex(i),:);
    matchingDataset{3,1}(i,:) = originalDataset{1,1}{3,1}(sampleIndex(i),:);
    matchingDataset{4,1}(i,:) = originalDataset{1,1}{4,1}(sampleIndex(i),:);
    matchingDataset{5,1}(i,:) = originalDataset{1,1}{5,1}(sampleIndex(i),:);
    matchingDataset{6,1}(i,:) = originalDataset{1,1}{6,1}(sampleIndex(i),:);
    
end

end

