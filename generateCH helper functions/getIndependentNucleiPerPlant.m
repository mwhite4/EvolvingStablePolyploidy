%% Martin White November 2018

%Helper function for generateComponentHomologs.  Decomposes data into that
%for individual nuclei of individual plants

%Ref: XXXXXXX

function [output] = getIndependentNucleiPerPlant(input)

[~,plants] = size(input);

for plant = 1:plants                                                        %For each independent plant
    
    [ID,~,IDindx]  = unique(input{6,plant});                                %Get the cell IDs
    IDcounts       = accumarray(IDindx,1);
    
    for cell = 1:length(ID)                                                 %For each unique cell ID,
        [~,LandCOs]     = size(input{1,plant});
        [~,SPSsites]    = size(input{2,plant});
        
        output{1,plant}{1,cell}     = nan(IDcounts(cell),LandCOs);
        output{1,plant}{2,cell}     = nan(IDcounts(cell),SPSsites);
        output{1,plant}{3,cell}     = nan(IDcounts(cell),LandCOs);          %Record of partner component homolog ID for each CO
        output{1,plant}{4,cell}     = strings(IDcounts(cell),SPSsites);     %Record of which homologs switched at each SPS site
        output{1,plant}{5,cell}     = strings(IDcounts(cell),1);            %Identifies the plant the homolog came from
        output{1,plant}{6,cell}     = nan(IDcounts(cell),1);
        
        for j = 1:length(input{6,plant})                                    %For every component homolog,
            if input{6,plant}(j,1) == ID(cell)                              %If it has the 'nth' cell ID,
                oIndex = sum(~isnan(output{1,plant}{1,cell}(:,1)))+1;
                output{1,plant}{1,cell}(oIndex,:) = input{1,plant}(j,:);    %put all of its data in column n of the output
                output{1,plant}{2,cell}(oIndex,:) = input{2,plant}(j,:);
                output{1,plant}{3,cell}(oIndex,:) = input{3,plant}(j,:);
                output{1,plant}{4,cell}(oIndex,:) = input{4,plant}(j,:);
                output{1,plant}{5,cell}(oIndex,:) = input{5,plant}(j,:);
                output{1,plant}{6,cell}(oIndex,:) = input{6,plant}(j,:);
            end
        end      
    end
    
end
end
