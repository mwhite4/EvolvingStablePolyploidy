function output = generateLengthBinnedCoC(input,lengthBins,numOfIntervals)

%% Martin White Jan 2020

% REF: XXXX

%Input: use the cell arrays created by generateComponentHomologs.m

%Output: a cell array.  
    %Each row is a different bin of chromosome lengths
    %column 1 is the minimum chromosome length for that bin
    %column 2 is the maxium chromosome length for that bin
    %column 1 is the number of chromosomes in that bin
    %column 3 is a cell array containing the CoC data for that bin, column
    %1 of which is inter-interval distance and column 2 is average
    %coefficient of coincidence

%%


%step 1: clean up input
input(input==0) = NaN;
input           = input(~all(isnan(input),2),:);
input           = input(:,~all(isnan(input)));

%Step 2: Bin measured chromosomes based on length
totalBins                   = length(lengthBins);
binnedChroms{totalBins-1,1} =[];
[totalChroms,~]             = size(input);
for chromosome = 1:totalChroms
    for bin = 1:totalBins-1
        if input(chromosome,1) >= lengthBins(bin) && ...
                input(chromosome,1) <= lengthBins(bin+1)
            binnedChroms{bin,1}(end+1,:) = input(chromosome,:);
        end
    end
end


%Step 3: generate CoC curves for each length bin
totalBins   = length(lengthBins);
output{totalBins-1,3} = [];

for bin = 1:totalBins-1
    output{bin,1} = lengthBins(bin);
    output{bin,2} = lengthBins(bin+1);
    output{bin,3} = length(binnedChroms{bin});
    if ~isempty(binnedChroms{bin})
        output{bin,4} = generateCoCcurves(binnedChroms{bin},numOfIntervals);
    end
end

    function output = generateCoCcurves(input,NumOfIntervals)
        
        %Step 1: Normalize CO positions relative to total chromosome length
        normCO = input(:,2:end)./input(:,1);
        [totChroms,~]   = size(normCO);
        
        
        %Step 2: Generate CoC data for each total number of intervals
        
        %Step 2.1: Calculate Interval Boundaries and labels
        intSize     = 1/NumOfIntervals;
        edges       = 0:intSize:1;
        
        %Step 2.2: Calculate the CO frequency for each interval
        obsCOfreq  = histcounts(normCO,edges);
        obsCOfreq  = obsCOfreq./totChroms;
        
        %Step 2.3: For each pair of intervals, calculate the observed
        %frequency of chromosomes with a(t least one) CO in both intervals
        COcounts(1:totChroms,1:NumOfIntervals) = nan;
        for chrom = 1:totChroms
            COcounts(chrom,:) = histcounts(normCO(chrom,:),edges);
        end
        max1COcounts    = COcounts;
        max1COcounts(max1COcounts>1) = 1;
        for j=1:NumOfIntervals-1
            for k=j+1:NumOfIntervals
                obs2COfreq(j,k-j)=(sum(max1COcounts(:,j).*...               %each row is an interval, each column is an inter-interval distance
                    max1COcounts(:,k)))./totChroms;
            end
        end
        
        %Step 2.4: For each pair of intervals, calculate the expected
        %frequency of double COs, assuming indepedence
        for j = 1:NumOfIntervals-1
            for k = j+1:NumOfIntervals
                exp2COfreq(j,k-j)  = obsCOfreq(j)*obsCOfreq(k);             %each row is an interval, each column is an inter-interval distance
            end
        end
        
        %Step 2.5: For each pair of intervals, calculate the Coefficient of
        %Coincidence
        CoC = obs2COfreq./exp2COfreq;
        
        %Step 2.6: Calculate average CoC for all inter-interval distances
        meanCoCperDistance   = nanmean(CoC,1);
        
        %Step 2.7: Calculate inter-interval distances
        interInterval           = 1:NumOfIntervals-1;
        relInterInterval        = 1/NumOfIntervals*interInterval;
        meanLength              = mean(input(:,1));
        interIntervalDist       = relInterInterval*meanLength;
        
        output(:,1) = interIntervalDist;
        output(:,2) = meanCoCperDistance;
        
    
    end

end

