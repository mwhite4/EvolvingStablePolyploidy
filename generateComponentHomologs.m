function [output] = generateComponentHomologs(lengthAndCOs,SPSsites,plantID,cellID)

%% Martin White November 2018

%This function takes the length, crossover (CO) positions, and synaptic
%parter switch (SPS) site positions along measured pachytene chromosome 
%bodies (bivalents and quadrivalents) and generates component homologs. It 
%also stores the plant and cell IDs of each chromosome for downstream
%processing.

%Helper Functions
    %This function requires the following six functions to run:
        %getEventPositions
        %getComponentHomologs
        %randomizeCHorientation
        %getIndependentPlants
        %getIndependentNucleiPerPlant
        %getCOCategories

%Ref: XXXXXXX

%Input

%NaNs are used to fill blank 'cells'

%lengthAndCOs: for each measured chromosome, the chromosome length is in
%column 1 and CO 'distances' are in subsequent columns.  'Distances' are
%the distance of the 'first CO' from the 'first' chromosome end, and then
%the distance of the next CO from the previous CO. Chromosome (genetic)
%identity are unknown.  All chromosomes within each nucleus was measured.

%SPSsites: contains SPS 'distances', measured analogously to CO distances
%(above), and from the same chromosome end.

%plantID: contains a string used to identify the source plant material

%cellID: contains an integer used to identify which cell of the plant the
%chromosome belonged to.


%Output
%Output takes the same format as the input, except that there are two
%additional variables and every row is a component homolog rather than a
%measured chromosome.

%outputFinal{3,1} contains the original homolog identity (1,2, 3, or 4)
%in column 1; outputFinal{3,1}(1:n,1).  Columns two onwards contains the
%identity of the pairing partner for each CO; outputFinal{3,1}(1:n,2:p)

%outputFinal{4,1} contains information on which of the two top and bottom
%homologs switched at each SPS site.

%outputFinal{5,1} identifies the plant the homolog came from (plant ID, a 
%string)

%outputFinal{6,1} identifies the nucleus the homolog came from (cell ID,
%integer)
%%

%Step 1: organize input data
lengthAndCOs(lengthAndCOs==0)   = NaN;                                      %convert any 0s to NaNs.  Blanks need to be NaNs,
SPSsites(SPSsites==0)           = NaN;                                      %however they tend to become 0s when imported from excel
input{1,1} = lengthAndCOs;
input{2,1} = SPSsites;
input{3,1} = plantID;
input{4,1} = cellID;

%Step 2: convert distance measurements to CO/SPS positions 
outputPositions = getEventPositions(input);

%Step 3: convert measured chromosomes to component homologs
componentHomologs(:,1) = getComponentHomologs(outputPositions);
componentHomologs(:,2) = getComponentHomologs(outputPositions);
componentHomologs(:,3) = getComponentHomologs(outputPositions);

%Step 4: re-orient component homologs, 3 independent times for
%component homolog dataset 1
reorientedCHs   = randomizeCHorientation(componentHomologs);
CH1_RO(:,1)     = reorientedCHs(:,1);
RO2             = randomizeCHorientation(componentHomologs);
RO3             = randomizeCHorientation(componentHomologs);

CH1_RO(:,2)     = RO2(:,1);
CH1_RO(:,3)     = RO3(:,1);


%Step 5: get CH1 RO1 data for each plant line
RO1CH1_Plant = getIndependentPlants(CH1_RO(:,1));

%Step 6: get CH1 RO1 data for each nucleus for each plant line
RO1CH1_PlantNuclei = getIndependentNucleiPerPlant(RO1CH1_Plant);

%Step 7: get E0, single CO, single partner multi-CO and multi-partner
%multi-CO categories
RO1CH1_COcategories = getCOCategories(CH1_RO(:,1));

%Step 8: format output
output{1,1} = CH1_RO;                                                       %first CH dataset, 3 different orientations
output{1,2} = reorientedCHs;                                                %3 different CH datasets (all reoriented same way)
output{1,3} = RO1CH1_Plant;                                                 %individual plants
output{1,4} = RO1CH1_PlantNuclei;                                           %individual nuclei per plant
output{1,5} = RO1CH1_COcategories;                                          %CO categories

end

