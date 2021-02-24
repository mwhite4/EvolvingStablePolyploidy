%% Martin White November 2018

%Helper function for generateComponentHomologs.  Flips or not, the 
%orientation of pachytene bivalents and quadrivalents with a 50%
%probability

%Ref: XXXXXXX


function [output] = randomizeCHorientation(input)

output = input;

[~,CHs]     = size(input);                                                  %Total number of independent component homolog datasets
[mCO,nCO]   = size(input{1,1});                                             %Total number of component homologs (should be same in each independent dataset), maximum number of COs
[~,nSPS]    = size(input{2,1});                                             %Maximum number of SPS sites

inputIndex  = 1;
outputIndex = 1;

while inputIndex <mCO+1                                                     %For every component homolog.
    
    if isnan(input{2,1}(inputIndex,1))                                      %If homolog is part of a synaptic bivalent
        x = rand(1);
        if x<0.5                                                            %flip this and its partner homolog with a 50% probability
            for dataset = 1:CHs
                output{1,dataset}(outputIndex,2:nCO)   = input{1,dataset}(inputIndex,1)   - input{1,dataset}(inputIndex,2:nCO);
                output{1,dataset}(outputIndex+1,2:nCO) = input{1,dataset}(inputIndex+1,1) - input{1,dataset}(inputIndex+1,2:nCO);
            end
        end
        outputIndex = outputIndex+2;
        inputIndex  = inputIndex+2;
        
    else                                                                    %if homolog is part of a synaptic quadrivalent
        x = rand(1);
        if x<0.5                                                            %flip this and its 3 partner homologs with a 50% probability
            for dataset = 1:CHs
                output{1,dataset}(outputIndex,2:nCO)   = input{1,dataset}(inputIndex,1)   - input{1,dataset}(inputIndex,2:nCO);
                output{2,dataset}(outputIndex,1:nSPS)  = input{1,dataset}(inputIndex,1)   - input{2,dataset}(inputIndex,1:nSPS);
                output{1,dataset}(outputIndex+1,2:nCO) = input{1,dataset}(inputIndex+1,1) - input{1,dataset}(inputIndex+1,2:nCO);
                output{2,dataset}(outputIndex+1,1:nSPS)= input{1,dataset}(inputIndex+1,1) - input{2,dataset}(inputIndex+1,1:nSPS);
                output{1,dataset}(outputIndex+2,2:nCO) = input{1,dataset}(inputIndex+2,1) - input{1,dataset}(inputIndex+2,2:nCO);
                output{2,dataset}(outputIndex+2,1:nSPS)= input{1,dataset}(inputIndex+2,1) - input{2,dataset}(inputIndex+2,1:nSPS);
                output{1,dataset}(outputIndex+3,2:nCO) = input{1,dataset}(inputIndex+3,1) - input{1,dataset}(inputIndex+3,2:nCO);
                output{2,dataset}(outputIndex+3,1:nSPS)= input{1,dataset}(inputIndex+3,1) - input{2,dataset}(inputIndex+3,1:nSPS);
            end
        end
        outputIndex = outputIndex+4;
        inputIndex  = inputIndex+4;
    end
    
end

for dataset = 1:CHs                                                         %for every independent component homolog dataset
    
    for i=1:mCO                                                             %for every component homolog in that dataset
        [output{1,dataset}(i,2:nCO),COindex]  = sort(output{1,dataset}(i,2:nCO));     %put the event positions in ascending order
        [output{2,dataset}(i,:),SPSindex]     = sort(output{2,dataset}(i,:));
        output{3,dataset}(i,2:nCO)            = output{3,dataset}(i,COindex+1);
        output{4,dataset}(i,:)                = output{4,dataset}(i,SPSindex);
    end
end

end