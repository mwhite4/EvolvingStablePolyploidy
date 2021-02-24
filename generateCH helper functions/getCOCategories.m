%% Martin White November 2018

%Helper function for generateComponentHomologs.  Decomposes data into the
%following CO categories: all CHs, CHs with 0 COs, CHs with 1 CO, CHs with
%2-homolog double COs only, all with the same partner homolog, and CHS with
%1 or more 3 homolog double CO

%Ref: XXXXXXX

function [output] = getCOCategories(input)
%input is a single reoriented component homolog dataset

[CHs,LandCOs]   = size(input{1,1});
[~,SPSsites]    = size(input{2,1});

for i = 1:5                                                                 %For each category
    output{1,i}     = nan(1,LandCOs);
    output{2,i}     = nan(1,SPSsites);
    output{3,i}     = nan(1,LandCOs);                                       %Record of partner component homolog ID for each CO
    output{4,i}     = strings(1,SPSsites);                                  %Record of which homologs switched at each SPS site
    output{5,i}     = strings(1,1);                                         %Identifies the plant the homolog came from
    output{6,i}     = nan(1,1);
end

output(:,1) = input(:,1);                                                   %Put data for all CHs in column 1

for i = 1:CHs                                                               %For every component homolog
    if sum(~isnan(input{1,1}(i,:))) == 1                                    %If component chromosome has no COs
        oIndex                  = sum(~isnan(output{1,2}(:,1)))+1;          %Put its data in column 2 of output
        output{1,2}(oIndex,:)   = input{1,1}(i,:);
        output{2,2}(oIndex,:)   = input{2,1}(i,:);
        output{3,2}(oIndex,:)   = input{3,1}(i,:);
        output{4,2}(oIndex,:)   = input{4,1}(i,:);
        output{5,2}(oIndex,:)   = input{5,1}(i,:);
        output{6,2}(oIndex,:)   = input{6,1}(i,:);
        
    elseif sum(~isnan(input{1,1}(i,:))) == 2                                %If component homolog has a single CO
        oIndex                  = sum(~isnan(output{1,3}(:,1)))+1;          %Put its data in column 3 of output
        output{1,3}(oIndex,:)   = input{1,1}(i,:);
        output{2,3}(oIndex,:)   = input{2,1}(i,:);
        output{3,3}(oIndex,:)   = input{3,1}(i,:);
        output{4,3}(oIndex,:)   = input{4,1}(i,:);
        output{5,3}(oIndex,:)   = input{5,1}(i,:);
        output{6,3}(oIndex,:)   = input{6,1}(i,:);
        
    else
        %check condition
        if length(unique(input{3,1}(i,2:sum(~isnan(input{3,1}(i,:))))))==1  %If more than one CO, and all CO are with the same partner compoment homolog
            oIndex                  = sum(~isnan(output{1,4}(:,1)))+1;      %Put its data in column 4 of output
            output{1,4}(oIndex,:)   = input{1,1}(i,:);
            output{2,4}(oIndex,:)   = input{2,1}(i,:);
            output{3,4}(oIndex,:)   = input{3,1}(i,:);
            output{4,4}(oIndex,:)   = input{4,1}(i,:);
            output{5,4}(oIndex,:)   = input{5,1}(i,:);
            output{6,4}(oIndex,:)   = input{6,1}(i,:);
            
        else
            oIndex                  = sum(~isnan(output{1,5}(:,1)))+1;      %If more than one CO, and not all CO are with the same partner compoment homolog
            output{1,5}(oIndex,:)   = input{1,1}(i,:);                      %Put its data in column 5 of output
            output{2,5}(oIndex,:)   = input{2,1}(i,:);
            output{3,5}(oIndex,:)   = input{3,1}(i,:);
            output{4,5}(oIndex,:)   = input{4,1}(i,:);
            output{5,5}(oIndex,:)   = input{5,1}(i,:);
            output{6,5}(oIndex,:)   = input{6,1}(i,:);        
        end
    end
end

