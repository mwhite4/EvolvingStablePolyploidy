%% getComponentChromosomes; Martin White November 2018
%This function takes the length, crossover (CO) positions, and synaptic
%parter switch (SPS) site positions along measured chromosomes and
%generates component chromosomes (2 per measured chromosome).

%Ref: XXXXXXX

%Input
%input is a cell array with 'CO positions' in row 1 column 1 {1,1}, SPS
%postions in row 2 column 1 {2,1}, plant line ID in row 3 column 1 {3,1},
%and cell ID in row 4 column 1 {4,1}
%within these cells, every row is a separate measured chromosome (n) 
%{1:2,1}(1:n,:)
%measured chromosome length is in column 1 of 'CO positions' {1,1}(1:n,1)
%NaNs are used to fill to empty 'cells' in each matrix.
%NB. plant and cell IDs are not required to generate component chromosomes,
%but are required for downstream analyses.

%Output
%Output takes the same format as the input, except that there are two
%additional variables and every row is a component chromosome rather than a
%measured chromosome.

%outputFinal{3,1} contains the original chromosome identity (1,2, 3, or 4)
%in column 1; outputFinal{3,1}(1:n,1).  Columns two onwards contains the
%identity of the pairing partner for each CO; outputFinal{3,1}(1:n,2:p)

%outputFinal{4,1} contains information on which of the two top and bottom
%chromosomes switched at each SPS site.

%outputFinal{5,1} identifies the plant the chromosome came from (plant line
% ID, a string)

%outputFinal{6,1} identifies the nucleus the chromosome came from (cell ID,
%integer)

%%
function [outputFinal] = getComponentHomologs(input)

[mCO,nCO]       = size(input{1,1});                                         %Generate output matrices
[~,nSPS]        = size(input{2,1});

output{1,1}     = nan(2*mCO,2*(nCO-1)+1);                                   %Container for Cell length and CO positions
output{2,1}     = nan(2*mCO,nSPS);                                          %Container for SPS positions
output{3,1}     = nan(2*mCO,2*(nCO-1)+1);                                   %Container for initial chromosome identifier and for each CO, the identity of its partner component chromosome
output{4,1}     = strings(2*mCO,nSPS);                                      %Container to record which chromosomes switched at each SPS site
output{5,1}     = strings(2*mCO,1);                                         %Identifies the plant the measured chromosome came from from downstream analyses
output{6,1}     = nan(2*mCO,1);                                             %Identifies the nucleus the measured chromosmoes came from from downstream analyses

iIndex          = 1;
oIndex          = 1;

while iIndex < mCO+1
    if sum(~isnan(input{2,1}(iIndex,:))) == 0                               %if chromosome is part of a synaptic bivalent (only length in SPS matrix)
        
        output{1,1}(oIndex,1:nCO)     = input{1,1}(iIndex,1:nCO);           %CC1 is the measured chromosome
        output{3,1}(oIndex,1)         = 1;
        output{3,1}(oIndex,2:sum(~isnan(input{1,1}(iIndex,:))))     = 2;
        output{5,1}(oIndex,1)         = strtrim(input{3,1}(iIndex,1));      %removes any leading and lagging spaces from plant line ID
        output{6,1}(oIndex,1)         = input{4,1}(iIndex,1);
        
        output{1,1}(oIndex+1,1:nCO)   = input{1,1}(iIndex,1:nCO);           %CC2 is a duplicate of CC1
        output{3,1}(oIndex+1,1)       = 2;
        output{3,1}(oIndex+1,2:sum(~isnan(input{1,1}(iIndex,:))))   = 1;
        output{5,1}(oIndex+1,1)       = strtrim(input{3,1}(iIndex,1));      %removes any leading and lagging spaces from plant line ID
        output{6,1}(oIndex+1,1)       = input{4,1}(iIndex,1);
        
        oIndex = oIndex+2;                                                  %move to the next measured chromosome
        iIndex  = iIndex+1;
        
    else                                                                    %if chromosome is part of a synaptic quadrivalent
        
        CC1co       = nan(1,2*(nCO-1)+1);                                   %These are temporary containers for values
        CC1sps      = nan(1,nSPS);                                          %CC1, 2, 3, 4 refers to the order of component chromosomes.
        CC1id       = nan(1,2*(nCO-1)+1);                                   %As chromosomes switch partners, their designation as CC1, 2, etc. changes,
        CC1id(1)    = 1;                                                    %but their initial ID (1, 2, 3 or 4) remains constant.
        %Switching always occurs between 1 of the 'top' chromosomes (CC1 or CC2)
        CC2co       = nan(1,2*(nCO-1)+1);                                   %and one of the 'bottom' chromosomes (CC3 or CC4)
        CC2sps      = nan(1,nSPS);
        CC2id       = nan(1,2*(nCO-1)+1);
        CC2id(1)    = 2;
        
        CC3co       = nan(1,2*(nCO-1)+1);
        CC3sps      = nan(1,nSPS);
        CC3id       = nan(1,2*(nCO-1)+1);
        CC3id(1)    = 3;
        
        CC4co       = nan(1,2*(nCO-1)+1);
        CC4sps      = nan(1,nSPS);
        CC4id       = nan(1,2*(nCO-1)+1);
        CC4id(1)    = 4;
        
        CC1co(1)  = input{2,1}(iIndex,1);                                   %assign first segment of 'top' chromosome to CC1 and CC2
        CC1sps(1) = input{2,1}(iIndex,1);
        CC2co(1)  = input{2,1}(iIndex,1);
        CC2sps(1) = input{2,1}(iIndex,1);
        for i = 2:nCO
            if input{1,1}(iIndex,i) < input{2,1}(iIndex,1)
                CC1co(i) = input{1,1}(iIndex,i);
                CC1id(i) = CC2id(1);
                CC2co(i) = input{1,1}(iIndex,i);
                CC2id(i) = CC1id(1);
            end
        end
        
        CC3co(1)  = input{2,1}(iIndex+1,1);                                 %assign first segment of 'bottom' chromosome to CC3 and CC4
        CC3sps(1) = input{2,1}(iIndex+1,1);
        CC4co(1)  = input{2,1}(iIndex+1,1);
        CC4sps(1) = input{2,1}(iIndex+1,1);
        for i = 2:nCO
            if input{1,1}(iIndex+1,i) < input{2,1}(iIndex+1,1)
                CC3co(i) = input{1,1}(iIndex+1,i);
                CC3id(i) = CC4id(1);
                CC4co(i) = input{1,1}(iIndex+1,i);
                CC4id(i) = CC3id(1);
            end
        end
        
        output{1,1}(oIndex,:)      = CC1co;                                 %Do your first switch
        output{1,1}(oIndex+1,:)    = CC3co;                                 %CC1 and CC2 are identical to each other.  As are CC3 and CC4.
        output{1,1}(oIndex+2,:)    = CC2co;                                 %Either CC1 or CC2 switches with either CC3 or CC4
        output{1,1}(oIndex+3,:)    = CC4co;                                 %The actual chromosome identities are unknown
        
        output{2,1}(oIndex,:)      = CC1sps;                                %There are 4 permutations for each switch:
        output{2,1}(oIndex+1,:)    = CC3sps;                                %CC1 and CC3, CC1 and CC4, CC2 and CC3, or CC3 and CC4 could switch
        output{2,1}(oIndex+2,:)    = CC2sps;                                %In this first case, since the actual identities are unknown,
        output{2,1}(oIndex+3,:)    = CC4sps;                                %all 4 permutations give the same result, so for ease, I just pick one:
        %switch CC2 and CC3 (its easier to draw in a diagram)
        output{3,1}(oIndex,:)      = CC1id;
        output{3,1}(oIndex+1,:)    = CC3id;
        output{3,1}(oIndex+2,:)    = CC2id;
        output{3,1}(oIndex+3,:)    = CC4id;
        
        output{4,1}(oIndex+1,1)    = {'switch CC2 & CC3'};
        output{4,1}(oIndex+2,1)    = {'switch CC2 & CC3'};
        
        
        if sum(~isnan(input{2,1}(iIndex,:))) >1                             %If there is more than 1 SPS
            for i = 2:sum(~isnan(input{2,1}(iIndex,:)))                     %For the second SPS onwards:
                
                CC1co  = output{1,1}(oIndex,:);                             %Re-define CCs after the previous switch events
                CC1sps = output{2,1}(oIndex,:);
                CC1id  = output{3,1}(oIndex,:);
                CC2co  = output{1,1}(oIndex+1,:);
                CC2sps = output{2,1}(oIndex+1,:);
                CC2id  = output{3,1}(oIndex+1,:);
                CC3co  = output{1,1}(oIndex+2,:);
                CC3sps = output{2,1}(oIndex+2,:);
                CC3id  = output{3,1}(oIndex+2,:);
                CC4co  = output{1,1}(oIndex+3,:);
                CC4sps = output{2,1}(oIndex+3,:);
                CC4id  = output{3,1}(oIndex+3,:);
                
                for j = 2:nCO                                               %Assign current segment of 'top' measured chromosome to current CC1 and CC2
                    if input{1,1}(iIndex,j) >input{2,1}(iIndex,i-1) && input{1,1}(iIndex,j) <input{2,1}(iIndex,i)
                        CC1co(sum(~isnan(CC1co(:)))+1) = CC1co(1) + input{1,1}(iIndex,j) - input{2,1}(iIndex,i-1);
                        CC1id(sum(~isnan(CC1id(:)))+1) = CC2id(1);
                        CC2co(sum(~isnan(CC2co(:)))+1) = CC2co(1) + input{1,1}(iIndex,j) - input{2,1}(iIndex,i-1);
                        CC2id(sum(~isnan(CC2id(:)))+1) = CC1id(1);
                    end
                end
                CC1sps(i) = CC1co(1) + input{2,1}(iIndex,i) - input{2,1}(iIndex,i-1);
                CC1co(1)  = CC1sps(i);
                CC2sps(i) = CC2co(1) + input{2,1}(iIndex,i) - input{2,1}(iIndex,i-1);
                CC2co(1)  = CC2sps(i);
                
                for j = 2:nCO                                               %Assign current segment of 'bottom' measured chromosome to current CC3 and CC4
                    if input{1,1}(iIndex+1,j) >input{2,1}(iIndex+1,i-1) && input{1,1}(iIndex+1,j) <input{2,1}(iIndex+1,i)
                        CC3co(sum(~isnan(CC3co(:)))+1) = CC3co(1) + input{1,1}(iIndex+1,j) - input{2,1}(iIndex+1,i-1);
                        CC3id(sum(~isnan(CC3id(:)))+1) = CC4id(1);
                        CC4co(sum(~isnan(CC4co(:)))+1) = CC4co(1) + input{1,1}(iIndex+1,j) - input{2,1}(iIndex+1,i-1);
                        CC4id(sum(~isnan(CC4id(:)))+1) = CC3id(1);
                    end
                end
                CC3sps(i) = CC3co(1) + input{2,1}(iIndex+1,i) - input{2,1}(iIndex+1,i-1);
                CC3co(1)  = CC3sps(i);
                CC4sps(i) = CC4co(1) + input{2,1}(iIndex+1,i) - input{2,1}(iIndex+1,i-1);
                CC4co(1)  = CC4sps(i);
                
                x = rand(1);                                                %Decide at random, which of the 4 possible switch events to do;
                y = rand(1);                                                %each with equal probability
                
                if x <0.5 && y <0.5                                         %switch CC1 and CC3
                    output{1,1}(oIndex,:)      = CC3co;
                    output{1,1}(oIndex+1,:)    = CC2co;
                    output{1,1}(oIndex+2,:)    = CC1co;
                    output{1,1}(oIndex+3,:)    = CC4co;
                    output{2,1}(oIndex,:)      = CC3sps;
                    output{2,1}(oIndex+1,:)    = CC2sps;
                    output{2,1}(oIndex+2,:)    = CC1sps;
                    output{2,1}(oIndex+3,:)    = CC4sps;
                    output{3,1}(oIndex,:)      = CC3id;
                    output{3,1}(oIndex+1,:)    = CC2id;
                    output{3,1}(oIndex+2,:)    = CC1id;
                    output{3,1}(oIndex+3,:)    = CC4id;
                    output{4,1}(oIndex,i)      = {'switch CC1 & CC3'};
                    output{4,1}(oIndex+2,i)    = {'switch CC1 & CC3'};
                    
                elseif x <0.5 && y >=0.5                                    %Switch CC1 and CC4
                    output{1,1}(oIndex,:)      = CC4co;
                    output{1,1}(oIndex+1,:)    = CC2co;
                    output{1,1}(oIndex+2,:)    = CC3co;
                    output{1,1}(oIndex+3,:)    = CC1co;
                    output{2,1}(oIndex,:)      = CC4sps;
                    output{2,1}(oIndex+1,:)    = CC2sps;
                    output{2,1}(oIndex+2,:)    = CC3sps;
                    output{2,1}(oIndex+3,:)    = CC1sps;
                    output{3,1}(oIndex,:)      = CC4id;
                    output{3,1}(oIndex+1,:)    = CC2id;
                    output{3,1}(oIndex+2,:)    = CC3id;
                    output{3,1}(oIndex+3,:)    = CC1id;
                    output{4,1}(oIndex,i)      = {'switch CC1 & CC4'};
                    output{4,1}(oIndex+3,i)    = {'switch CC1 & CC4'};
                    
                elseif x >=0.5 && y <0.5                                    %Switch CC2 and CC3
                    output{1,1}(oIndex,:)      = CC1co;
                    output{1,1}(oIndex+1,:)    = CC3co;
                    output{1,1}(oIndex+2,:)    = CC2co;
                    output{1,1}(oIndex+3,:)    = CC4co;
                    output{2,1}(oIndex,:)      = CC1sps;
                    output{2,1}(oIndex+1,:)    = CC3sps;
                    output{2,1}(oIndex+2,:)    = CC2sps;
                    output{2,1}(oIndex+3,:)    = CC4sps;
                    output{3,1}(oIndex,:)      = CC1id;
                    output{3,1}(oIndex+1,:)    = CC3id;
                    output{3,1}(oIndex+2,:)    = CC2id;
                    output{3,1}(oIndex+3,:)    = CC4id;
                    output{4,1}(oIndex+1,i)    = {'switch CC2 & CC3'};
                    output{4,1}(oIndex+2,i)    = {'switch CC2 & CC3'};
                    
                else                                                        %Switch CC2 and CC4
                    output{1,1}(oIndex,:)      = CC1co;
                    output{1,1}(oIndex+1,:)    = CC4co;
                    output{1,1}(oIndex+2,:)    = CC3co;
                    output{1,1}(oIndex+3,:)    = CC2co;
                    output{2,1}(oIndex,:)      = CC1sps;
                    output{2,1}(oIndex+1,:)    = CC4sps;
                    output{2,1}(oIndex+2,:)    = CC3sps;
                    output{2,1}(oIndex+3,:)    = CC2sps;
                    output{3,1}(oIndex,:)      = CC1id;
                    output{3,1}(oIndex+1,:)    = CC4id;
                    output{3,1}(oIndex+2,:)    = CC3id;
                    output{3,1}(oIndex+3,:)    = CC2id;
                    output{4,1}(oIndex+1,i)    = {'switch CC2 & CC4'};
                    output{4,1}(oIndex+3,i)    = {'switch CC2 & CC4'};
                end
            end
        end
        
        for j = 2:nCO                                                       %add final segment of 'top' measured chromosome to CC1 and CC2
            if input{1,1}(iIndex,j) >input{2,1}(iIndex,sum(~isnan(input{2,1}(iIndex,:))))
                output{1,1}(oIndex,sum(~isnan(output{1,1}(oIndex,:)))+1)     = output{1,1}(oIndex,1)   + input{1,1}(iIndex,j) - input{2,1}(iIndex,sum(~isnan(input{2,1}(iIndex,:))));
                output{3,1}(oIndex,sum(~isnan(output{3,1}(oIndex,:)))+1)     = output{3,1}(oIndex+1,1);
                output{1,1}(oIndex+1,sum(~isnan(output{1,1}(oIndex+1,:)))+1) = output{1,1}(oIndex+1,1) + input{1,1}(iIndex,j) - input{2,1}(iIndex,sum(~isnan(input{2,1}(iIndex,:))));
                output{3,1}(oIndex+1,sum(~isnan(output{3,1}(oIndex+1,:)))+1) = output{3,1}(oIndex,1);
            end
        end
        output{1,1}(oIndex,1)   = output{1,1}(oIndex,1)   + input{1,1}(iIndex,1) - input{2,1}(iIndex,sum(~isnan(input{2,1}(iIndex,:))));
        output{1,1}(oIndex+1,1) = output{1,1}(oIndex+1,1) + input{1,1}(iIndex,1) - input{2,1}(iIndex,sum(~isnan(input{2,1}(iIndex,:))));
        
        
        for j = 2:nCO                                                       %add final segment of 'bottom' measured chromosome to CC3 and CC4
            if input{1,1}(iIndex+1,j) >input{2,1}(iIndex+1,sum(~isnan(input{2,1}(iIndex+1,:))))
                output{1,1}(oIndex+2,sum(~isnan(output{1,1}(oIndex+2,:)))+1) = output{1,1}(oIndex+2,1) + input{1,1}(iIndex+1,j) - input{2,1}(iIndex+1,sum(~isnan(input{2,1}(iIndex+1,:))));
                output{3,1}(oIndex+2,sum(~isnan(output{3,1}(oIndex+2,:)))+1) = output{3,1}(oIndex+3,1);
                output{1,1}(oIndex+3,sum(~isnan(output{1,1}(oIndex+3,:)))+1) = output{1,1}(oIndex+3,1) + input{1,1}(iIndex+1,j) - input{2,1}(iIndex+1,sum(~isnan(input{2,1}(iIndex+1,:))));
                output{3,1}(oIndex+3,sum(~isnan(output{3,1}(oIndex+3,:)))+1) = output{3,1}(oIndex+2,1);
            end
        end
        output{1,1}(oIndex+2,1) = output{1,1}(oIndex+2,1) + input{1,1}(iIndex+1,1) - input{2,1}(iIndex+1,sum(~isnan(input{2,1}(iIndex+1,:))));
        output{1,1}(oIndex+3,1) = output{1,1}(oIndex+3,1) + input{1,1}(iIndex+1,1) - input{2,1}(iIndex+1,sum(~isnan(input{2,1}(iIndex+1,:))));
        
        output{5,1}(oIndex,1)           = strtrim(input{3,1}(iIndex,1));    %add the plant ID
        output{6,1}(oIndex,1)           = input{4,1}(iIndex,1);             %add the cell ID
        output{5,1}(oIndex+1,1)         = strtrim(input{3,1}(iIndex,1));    %NB. all component chromosomes have the same
        output{6,1}(oIndex+1,1)         = input{4,1}(iIndex,1);             %plant and cell IDs
        output{5,1}(oIndex+2,1)         = strtrim(input{3,1}(iIndex,1));
        output{6,1}(oIndex+2,1)         = input{4,1}(iIndex,1);
        output{5,1}(oIndex+3,1)         = strtrim(input{3,1}(iIndex,1));
        output{6,1}(oIndex+3,1)         = input{4,1}(iIndex,1);
        
        oIndex  = oIndex+4;
        iIndex  = iIndex+2;                                                 %skip to the next 'top'/pachytene bivalent measured chromosome
        
    end                                                                     

end                                                                         %get rid of excess NaNs in output matrices
    outputFinal{1,1} = output{1,1}(1:2*mCO,1:max(sum(~isnan(output{1,1}),2)));
    if max(sum(~isnan(output{2,1}),2))>0
        outputFinal{2,1} = output{2,1}(1:2*mCO,1:max(sum(~isnan(output{2,1}),2)));
    else
        outputFinal{2,1} = output{2,1};
    end
    outputFinal{3,1} = output{3,1}(1:2*mCO,1:max(sum(~isnan(output{3,1}),2)));
    outputFinal{4,1} = output{4,1};
    outputFinal{5,1} = output{5,1};
    outputFinal{6,1} = output{6,1};

end