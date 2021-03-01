
%% Martin White Jan 2020

% REF: XXXX

%function to generate modified CoC curves to compare single v
%multi-partner component homologs.

%Input: use the cell arrays created by generateComponentHomologs.m

%Output is a cell array.  

%Each row is a different bin of compononent homolog lengths.  

%Column 1 is minimum component homolog length
%Column 2 is maximum component homolog length
%Column 3 contains inter-interval distances and corresponding Coefficient
%of coincidence values for all component homologs
%Column 4 contains inter-interval distances and corresponding Coefficient
%of coincidence values for 2-homolog double crossover component homologs
%Column 5 contains inter-interval distances and corresponding Coefficient
%of coincidence values for 3-homolog double crossover component homologs

function [output] = generate2H3HLengthBinnedCoC(input,lengthBinEdges,numOfIntervals)


%Step 1: create containers for outputs
inputData       = input{1,5}(1,:);
totalBins       = length(lengthBinEdges) - 1;

obs2COfreq_2chrom{1,totalBins}          = [];
obs2COfreq_multiChrom{1,totalBins}      = [];
exp2COfreq{1,totalBins}                 = [];
CoC{1,totalBins}                        = [];
twoChromDoubleCoC{1,totalBins}          = [];
multiChromDoubleCoC{1,totalBins}        = [];
twoMeanCoCperDistance{1,totalBins}      = [];
multiMeanCoCperDistance{1,totalBins}    = [];
allMeanCoCperDistance{1,totalBins}      = [];
interIntervalDist{1,totalBins}          = [];
output{totalBins,5}                     = [];


%Step 2: Bin measured chromosomes based on length
binnedChroms{totalBins,5} =[];
for i = 1:5
    [totalChroms,~]             = size(inputData{1,i});   
    for chromosome = 1:totalChroms
        for bin = 1:totalBins
            if inputData{1,i}(chromosome,1) >= lengthBinEdges(bin) && ...
                    inputData{1,i}(chromosome,1) <= lengthBinEdges(bin+1)
                binnedChroms{bin,i}(end+1,:) = inputData{1,i}(chromosome,:);
            end
        end
    end
end


%For each length category:
for i = 1:totalBins
    
    %step 3.1:adjust the number of the 2-chromosome and multi-chromosome
    %double categories so that their freq match the freq of multi CO CCs in
    %the measured population
    
    allCCs                  = length(binnedChroms{i,1}(:,1));
    twoChromMultiCOCCs      = length(binnedChroms{i,4}(:,1));
    multiChromMultiCOCCs    = length(binnedChroms{i,5}(:,1));
    totalMultiCOCCs         = twoChromMultiCOCCs+multiChromMultiCOCCs;
    meanCClength            = mean(binnedChroms{i,1}(:,1));
    
    twoChromAdjPopn = round(twoChromMultiCOCCs/totalMultiCOCCs*allCCs);
    multiChromAdjPopn = round(multiChromMultiCOCCs/totalMultiCOCCs*allCCs);
    
    binnedChroms{i,4}(twoChromMultiCOCCs+1:twoChromAdjPopn,1) = meanCClength;
    binnedChroms{i,4}(twoChromMultiCOCCs+1:twoChromAdjPopn,2:end) = NaN;
    
    binnedChroms{i,5}(multiChromMultiCOCCs+1:multiChromAdjPopn,1) = meanCClength;
    binnedChroms{i,5}(multiChromMultiCOCCs+1:multiChromAdjPopn,2:end) = NaN; 
    
    
    %step 3.1: get the observed 2CO freq for all 2-chromosome doubles
    [obs2COfreq_2chrom(:,i),~,~] = getCoC(binnedChroms(i,4),numOfIntervals);
    
    %step 3.2: get the observed 2 CO freq for all multi-chromosome doubles
    [obs2COfreq_multiChrom(:,i),~,~] = getCoC(binnedChroms(i,5),numOfIntervals);
    
    %step 3.3: calculate the expected double CO freq using all CO cats
    [~,exp2COfreq(:,i),CoC(:,i)] = getCoC(binnedChroms(i,1),numOfIntervals);
    
    
    %Step 3.4: calculate the CoC values for 2-homolog double COs
    twoChromDoubleCoC{1,i} = obs2COfreq_2chrom{1,i}./exp2COfreq{1,i};
    %Step 3.5: calclulate the CoC values for 3-homolog double COs
    multiChromDoubleCoC{1,i} = obs2COfreq_multiChrom{1,i}./exp2COfreq{1,i};
    %Step 3.6: calculate the mean CoC values for 2-homolog double COs
    twoMeanCoCperDistance{1,i}   = nanmean(twoChromDoubleCoC{1,i},1);
    %Step 3.7: calclulate the CoC values for 3-homolog double COs
    multiMeanCoCperDistance{1,i} = nanmean(multiChromDoubleCoC{1,i},1);
    %Step 3.8: calclulate the CoC values for the whole population
    allMeanCoCperDistance{1,i}   = nanmean(CoC{1,i},1);
    %Step 3.9: calculate the inter-interval distances
    interInterval           = 1:numOfIntervals-1;
    relInterInterval        = 1/numOfIntervals*interInterval;
    interIntervalDist{1,i}  = relInterInterval*meanCClength;
    
    output{i,1}(:,1) = lengthBinEdges(i);
    output{i,2}(:,1) = lengthBinEdges(i+1);
    output{i,3}(:,1) = interIntervalDist{1,i}';
    output{i,3}(:,2) = allMeanCoCperDistance{1,i}';
    output{i,4}(:,1) = interIntervalDist{1,i}';
    output{i,4}(:,2) = twoMeanCoCperDistance{1,i}';
    output{i,5}(:,1) = interIntervalDist{1,i}';
    output{i,5}(:,2) = multiMeanCoCperDistance{1,i}';
    
end

%accessory function
function [obs2COfreq,exp2COfreq,CoC] = getCoC(input,NumOfIntervals)

normCO{1,length(input)}                 = [];
COcounts{1,length(input)}               = [];
obsCOfreq{1,length(input)}              = [];
obs2COfreq{1,length(input)}             = [];
exp2COfreq{1,length(input)}             = [];
CoC{1,length(input)}                    = [];
meanCoCperDistance{1,length(input)}     = [];
medianCoCperDistance{1,length(input)}   = [];
interIntervalDistance{1,length(input)}  = [];

for dataset = 1:length(input)
    
    %Step 1: Normalize CO positions relative to total chromosome length
    normCO{dataset} = input{dataset}(:,2:end)./input{dataset}(:,1);
    [totChroms,~]   = size(normCO{dataset});
    
    
    %Step 2: Generate CoC data for each total number of intervals
        
        %Step 2.1: Calculate Interval Boundaries and labels
        intSize     = 1/NumOfIntervals;
        edges       = 0:intSize:1;
               
        %Step 2.2: Calculate the CO frequency for each interval
        obsCOfreq{1,dataset}  = histcounts(normCO{dataset},edges);
        obsCOfreq{1,dataset}  = obsCOfreq{1,dataset}./totChroms;
        
        %Step 2.3: For each pair of intervals, calculate the observed
        %frequency of chromosomes with a(t least one) CO in both intervals
        COcounts{1,dataset}(1:totChroms,1:NumOfIntervals) = nan;
        for chrom = 1:totChroms
            COcounts{1,dataset}(chrom,:) = histcounts(normCO{dataset}(chrom,:),edges);
        end
        max1COcounts    = COcounts{1,dataset};
        max1COcounts(max1COcounts>1) = 1;
        for j=1:NumOfIntervals-1
            for k=j+1:NumOfIntervals
                obs2COfreq{1,dataset}(j,k-j)=(sum(max1COcounts(:,j).*...    %each row is an interval, each column is an inter-interval distance
                    max1COcounts(:,k)))./totChroms;
            end
        end
        
        %Step 2.4: For each pair of intervals, calculate the expected
        %frequency of double COs, assuming indepedence
        for j = 1:NumOfIntervals-1
            for k = j+1:NumOfIntervals(1)
                exp2COfreq{1,dataset}(j,k-j)  = obsCOfreq{1,dataset}(j)...  %each row is an interval, each column is an inter-interval distance
                    *obsCOfreq{1,dataset}(k);
            end
        end
        
        %Step 2.5: For each pair of intervals, calculate the Coefficient of
        %Coincidence
        CoC{1,dataset} = obs2COfreq{1,dataset}./exp2COfreq{1,dataset};
        
        %Step 2.6: Calculate average CoC for all inter-interval distances
        meanCoCperDistance{1,dataset}   = nanmean(CoC{1,dataset},1);
        medianCoCperDistance{1,dataset} = nanmedian(CoC{1,dataset},1);
        
        %Step 2.7: Calculate inter-interval distances
        inter_Interval      = 1:NumOfIntervals-1;
        relInter_Interval   = 1/NumOfIntervals*inter_Interval;
        meanLength          = mean(input{dataset}(:,1));
        interIntervalDistance{1,dataset} = relInter_Interval*meanLength;

end

end

end  
