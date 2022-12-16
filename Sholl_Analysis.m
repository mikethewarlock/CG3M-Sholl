%% CG3M Sholl analysis
% Mike Giannetto
% Nedergaard Lab
% 12/7/2022

%{
This code takes analysis done by PT (see his notebook) in FIJI with
thresholding, skeletonizing, and sholl measurements. Images were taken in
hippocampus, confocal imaging.

>>>>>>>>>>>>>>>INPUTS>>>>>>>>>>>>>>>>
excelFile = uigetfile; 
    This is an excel file where the column is distance
    from the center of the astrocyte (in micrometers), then each subsequent
    column is the number of intersections the skeletonaized astrocyte has
    with that ring. The first row is animal ID, the rest of the rows are
    data.
    If an animal does not have as many astrocytes as another, those columns
    are left as blank cells. Also, if an astrocyte does not have any
    intersections past a certain distance, those cells are also blank.

<<<<<<<<<<<<<<<OUTPUTS<<<<<<<<<<<<<<<
Graphs / excel file / data saved as desired (TBD as of 12/7/22)


*******************NOTES*************************
Using estimates from Oberheim paper: https://pubmed.ncbi.nlm.nih.gov/19279265/ 
Protoplasmic astrocytes ~ 60um diameter, with 4-6 main processes 
                            (more cortical)
Fibrous astrocytes ~ 80um diameter
                            (more in W matter, hippocampal)

Script is designed to work with 20 mice, will need to update any instances
of 20 with a variable "numAnimals" if using with other data, or manually
replace the number 20 with however many animals are in spreadsheet

NEED TO CONFIRM SCALE / STEP SIZE of each sholl ring in analysis. Will
likely need to update excel file and "scale" variable in the code

Update 12/16/22
Code now calculates mean across all astrocytes, and sum across all
astrocytes.

need to further analyze individual cells, ie how many processes. Sum total
looks quite different from mean, so there is some bimodal / weird effect
going on
%}

%% Code starts here: Load in excel file, organize data into cell array

% Scale here refers to size of steps between concentric sholl rings
% (radius of sholl ring)
stepSize = 2; %Steps in px used in sholl analysis
scale = (0.09655)* stepSize; %Scale is correct .09655 px/um 
numAnimals = str2double(inputdlg('How many total animals?'));

% dim1 IND for 100um is 324 (cell array) or 323(data)
% Scale - every row is .310742um

excelFile = uigetfile('*.xlsx'); %Get excel filename
[~, ~, cellData] = xlsread(excelFile); %Load file

%Create table of animal IDs
animalID = cell(numAnimals,1); 
for ii=1:numAnimals
    tmpNum = char(num2str(ii));
     tmpChar = ['Animal ' tmpNum];
     animalID{ii} = string(tmpChar);
end

%Convert column titles (first row in cellData) to string array
for ii = 1:size(cellData,2)
    cellData{1,ii} = string(cellData{1,ii});
end

%Create matrix of values, create distance value
numData = cell2mat(cellData(2:end,:));
dist = [1:size(numData,1)];
dist = dist .* scale ;
% Create a cell array, each cell is a matrix containing all the data for that animal
% Each animal can have varying size matrix (number of astrocytes)
shollRaw = cell(numAnimals,1); 

for ii = 1:numAnimals
    tmpMat = [];
    tmpTF = zeros(size(cellData,2),1);
    for kk = 1:size(cellData,2)
        tmpTF(kk) = strcmp(animalID{ii}, cellData{1,kk});
    end
    tmpIND = find(tmpTF); %find columns that match animal
    tmpMat = numData(:,tmpIND); %Create only number matrix
    tmpInd2= find(~isnan(tmpMat(2,:)));
    tmpMat = tmpMat(:,tmpInd2); %Remove NaN columns (no astrocyte)
    shollRaw{ii} = tmpMat;
end
%Convert NaN values in each column to 0
shollNanZ = cell(numAnimals,1);
for ii = 1:numAnimals
    tmpMat = shollRaw{ii};
    tmpCellNum = size(tmpMat,2);
    noNanMat = zeros(size(tmpMat));
    for kk = 1:tmpCellNum
        tmpVecNan = find(isnan(tmpMat(:,kk)));
        tmpVecZero = tmpMat(:,kk); tmpVecZero(tmpVecNan) = 0;
        noNanMat(:,kk) = tmpVecZero;
    end
    shollNanZ{ii} = noNanMat;
end


%Load in animal info from excel file "CG3M DEMOGRAPHICS"
% Should edit filename if using different data, but keep format consistent
[~,~,demo] = xlsread('CG3M demographics.xlsx');
grpInd = cell2mat(demo(2:(numAnimals+1),6));

%Save variables as a -mat fil "shollRaw.mat"
save ("shollRaw",'shollRaw','shollNanZ','dist','demo','grpInd','scale','numAnimals')

%% Quantifying Sholl Raw
clear
load 'shollRaw.mat'
% Determine furthest process with an intersection for each astrocyte
diameterAstro = cell(numAnimals,1);

for ii = 1:numAnimals
    tmpMat = shollNanZ{ii};
    tmpCellNum = size(tmpMat,2);
    tmpDiam = zeros(tmpCellNum,0);
    for kk = 1:tmpCellNum
        tmpDiam(kk) = find(tmpMat(:,kk),1,'last');
    end
    diameterAstro{ii} = tmpDiam;
end
clear tmpMat tmpCellNum tmpDiam 

    %Reorganize diameter data into single matrix, each column is an
    %astrocyte, order of astrocytes is animal 1, then animal 2 etc. Varying
    %number of astrocytes per animal means that will need to specify which
    %column belongs to which animal (Line 142, Line 155)
diamTogether=[];
for ii=1:numAnimals
    diamTogether = [diamTogether diameterAstro{ii}];
end
diamTogetherScale = diamTogether' .* scale; %scale index to um value

    %Count number of Astrocytes per Animal
numAstros = [];
for ii = 1:numAnimals
    numAstros(ii) = length(diameterAstro{ii});
end
numAstros = numAstros';

    %Create indexing for each animal when data is arranged as a single matrix
anmlTogInd = cell(numAnimals,1);
anmlTogInd{1} = (1:numAstros(1));
for ii = 2:numAnimals
    sumTmp1 = sum(numAstros(1:ii-1),'all')+1;
    sumTmp2 = sum(numAstros(1:ii),'all');
    anmlTogInd{ii} = [sumTmp1:sumTmp2];
end

    %Create variable for group ID to be used with cell array of 20 mice
grpAnimalIndSaline = find((grpInd .* -1)+1);
grpAnimalIndGlucose = find(grpInd);

numSal = length(grpAnimalIndSaline);
numGlu = length(grpAnimalIndGlucose);
    %Create variables for group ID, to be used with matrix containing all data
grpTogSal = [];
grpTogGlu = [];
for ii = 1:numSal
    grpTogSal = [grpTogSal anmlTogInd{grpAnimalIndSaline(ii)}];
end
for ii = 1:numGlu
    grpTogGlu = [grpTogGlu anmlTogInd{grpAnimalIndGlucose(ii)}];
end


figure %Plot diameter all data together
plot(diamTogetherScale,'.k')

figure %Plot saline vs glucose, animals pooled together
subplot(1,2,1)
plot(diamTogetherScale(grpTogSal),'.b')
ylim([0 100])
subplot(1,2,2)
plot(diamTogetherScale(grpTogGlu),'.r')
ylim([0 100])

figure %Plot individual animals, saline on top, glucose on bottom
for ii = 1:numSal
    subplot(2,numSal,ii)
    plot(diamTogetherScale(anmlTogInd{grpAnimalIndSaline(ii)}),'.b')
    ylim([0 200])
end
for ii= 1:numGlu
    subplot(2,numGlu,ii+numSal)
    plot(diamTogetherScale(anmlTogInd{grpAnimalIndGlucose(ii)}),'.r')
    ylim([0 200])
end

%% Make a histogram for each group based on furthest distance
binsHist1 = (0:5:100); % need to adjust if change scale, binning in um 

histDiamSal = histc(diamTogetherScale(grpTogSal),binsHist1);
histDiamGlu = histc(diamTogetherScale(grpTogGlu),binsHist1);

figure
subplot(2,1,1)
bar(histDiamSal,'b')
%xticks([0:2:20]);
%xticklabels(( [0:2:20] .*10))
ylim([0 25])
ylabel('Number of Astrocytes')
xlabel('Furthest Process(um)')
xticks([1:4:21]);
xticklabels(binsHist1(1:4:21));

subplot(2,1,2)
bar(histDiamGlu,'r')
%xticks([0:2:20]);
%xticklabels(( [0:2:20] .*10))
ylim([0 25])   
ylabel('Number of Astrocytes')
xlabel('Furthest Process(um)')
xticks([1:4:21]);
xticklabels(binsHist1(1:4:21));
%% Plot sholl curves (number of intersections)
%{ 
Going to use "shollRaw" (where values with 0 -> NaN), this will allow
exclusion of zero values when quantifying number of intersections, and
avoid having the averages dragged down by astrocytes with no processing
in outer bounds
%}
allAstros = [];
for ii = 1:numAnimals
    tmpMat = shollRaw{ii};
    allAstros = [allAstros tmpMat];
end

meanIntersectSal = mean(allAstros(:,grpTogSal),2,'omitnan');
stdevIntSal = std(allAstros(:,grpTogSal),0,2,'omitnan');
sampleSizeSal = sum(~isnan(allAstros(:,grpTogSal)), 2);
stdErrorSal = stdevIntSal ./ (sqrt(sampleSizeSal));

meanIntersectGlu = mean(allAstros(:,grpTogGlu),2,'omitnan');
stdevIntGlu = std(allAstros(:,grpTogGlu),0,2,'omitnan');
sampleSizeGlu = sum(~isnan(allAstros(:,grpTogGlu)),2);
stdErrorGlu = stdevIntGlu ./ (sqrt(sampleSizeGlu));

figure
plot(dist,meanIntersectSal);hold on;
plot(dist,meanIntersectGlu); hold off;
title('Mean intersections')

sumIntSal = sum(allAstros(:,grpTogSal),2,'omitnan');
sumIntGlu = sum(allAstros(:,grpTogGlu),2,'omitnan');
figure
plot(dist,sumIntSal);hold on;
plot(dist,sumIntGlu);hold off;
title('Sum Intersections')

sumIntSalIndiv = sum(allAstros(:,grpTogSal),1,'omitnan');
sumIntGluIndiv = sum(allAstros(:,grpTogGlu),1,'omitnan');
figure
plot(zeros(length(grpTogSal),1),sumIntSalIndiv,'.b',...
    'MarkerSize',20);hold on;
plot(ones(length(grpTogGlu),1),sumIntGluIndiv,'.r',...
    'MarkerSize',20);hold off;


% Need to figure out good way to make ghost error bars

% figure
% figHand = plot(dist,meanIntersectSal);
% tmpErrorBar = errorbar(dist,meanIntersectSal,stdErrorSal);
% tmpErrorBar

    
    
    