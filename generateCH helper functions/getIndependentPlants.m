%% Martin White November 2018

%Helper function for generateComponentHomologs.  Decomposes data into that
%for individual plants

%Ref: XXXXXXX

function [output] = getIndependentPlants(input)

[ID,~,IDindx]  = unique(input{5,1});                                       %Get the IDs of the different plants
IDcounts        = accumarray(IDindx,1);

for i = 1:length(ID)                                                       %For each unique plant ID,
    [~,LandCOs]     = size(input{1,1});
    [~,SPSsites]    = size(input{2,1});
    
    output{1,i}     = nan(IDcounts(i),LandCOs);
    output{2,i}     = nan(IDcounts(i),SPSsites);
    output{3,i}     = nan(IDcounts(i),LandCOs);                             %Record of partner component homolog ID for each CO
    output{4,i}     = strings(IDcounts(i),SPSsites);                        %Record of which homologs switched at each SPS site
    output{5,i}     = strings(IDcounts(i),1);                               %Identifies the plant the measured chromosome came from
    output{6,i}     = nan(IDcounts(i),1);
    
    for j = 1:length(input{5,1})                                            %For every component homolog,
        if input{5,1}(j,1) == ID(i)                                         %If it has the 'nth' ID,
            oIndex                = sum(~isnan(output{1,i}(:,1)))+1;
            output{1,i}(oIndex,:) = input{1,1}(j,:);                        %put all of its data in column n of the output
            output{2,i}(oIndex,:) = input{2,1}(j,:);
            output{3,i}(oIndex,:) = input{3,1}(j,:);
            output{4,i}(oIndex,:) = input{4,1}(j,:);
            output{5,i}(oIndex,:) = input{5,1}(j,:);
            output{6,i}(oIndex,:) = input{6,1}(j,:);                   
        end
    end
    
end
end