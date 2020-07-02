function obj = results2library(obj)

% Open report file
reportData = readReport(obj);

% Filter PSM
scores = reportData(2:end,24);
numericScores = [];
for j = 1:length(scores)
  numericScores = [numericScores;str2num(char(scores(j,1)))];
end
removeIDX = find(numericScores(:,1)<obj.settings.minPSMScore | isnan(numericScores(:,1)));
removeIDX = removeIDX+1;
reportData(removeIDX,:) = [];

% Get files and spectrum titles
sequenceList = reportData(2:end,3);
proteinList = reportData(2:end,2);
mgfList = unique(reportData(2:end,10));
scanList = reportData(2:end,11);

% Collect peptide spectra and store as .mat
cd(obj.folder.identification);
count = 0;
for j = 1:length(mgfList)
    MGFStruct = readMGF(mgfList{j});
    titleList = {MGFStruct.scan.scanName}';
    for n = 1:length(scanList)
        scanIndex = find(strcmp(titleList,scanList{n}));
        if ~isempty(scanIndex)
           count = count+1; 
           library(count).sequence = sequenceList{n};
           library(count).protein = proteinList{n};
           filteredSpectrum = retrieveAnnotations(sequenceList{n},MGFStruct.scan(scanIndex).scanData);
           library(count).spectrum = filteredSpectrum;
        end
    end
end
obj.output.library = library;
cd([obj.folder.identification '\library']);
save('library.mat','library');
cd(obj.folder.mainFolder);

end

function filteredSpectrum = retrieveAnnotations(sequence,scanData)
    [yseries,bseries] = fragmentSequence(sequence);
    annotatedIons = [];
    tolerance = 20; % ppm
    
    yions = yseries.mz;
    for j = 1:length(yions)
        maxDev = ppmDeviation(yions(j),tolerance);
        idx = find(scanData(:,1) > yions(j)-maxDev & scanData(:,1) < yions(j)+maxDev);
        if ~isempty(idx)
           if numel(idx) > 1
               diff = abs(scanData(idx,1)-yions(j));
               min_diff = find(diff==min(diff));
               annotatedIons = [annotatedIons;idx(min_diff)];
           else
               annotatedIons = [annotatedIons;idx];
           end
        else
           continue 
        end
    end
    
    bions = bseries.mz;
    for j = 1:length(bions)
        maxDev = ppmDeviation(bions(j),tolerance);
        idx = find(scanData(:,1) > bions(j)-maxDev & scanData(:,1) < bions(j)+maxDev);
        if ~isempty(idx)
           if numel(idx) > 1
               diff = abs(scanData(idx,1)-bions(j));
               min_diff = find(diff==min(diff));
               annotatedIons = [annotatedIons;idx(min_diff)];
           else
               annotatedIons = [annotatedIons;idx];
           end
        else
           continue 
        end
    end
    
    yWater = yseries.waterLoss;
    for j = 1:length(yWater)
        maxDev = ppmDeviation(yWater(j),tolerance);
        idx = find(scanData(:,1) > yWater(j)-maxDev & scanData(:,1) < yWater(j)+maxDev);
        if ~isempty(idx)
           if numel(idx) > 1
               diff = abs(scanData(idx,1)-yWater(j));
               min_diff = find(diff==min(diff));
               annotatedIons = [annotatedIons;idx(min_diff)];
           else
               annotatedIons = [annotatedIons;idx];
           end
        else
           continue 
        end
    end
    
    yAmmonia = yseries.ammoniaLoss;
    for j = 1:length(yAmmonia)
        maxDev = ppmDeviation(yAmmonia(j),tolerance);
        idx = find(scanData(:,1) > yAmmonia(j)-maxDev & scanData(:,1) < yAmmonia(j)+maxDev);
        if ~isempty(idx)
           if numel(idx) > 1
               diff = abs(scanData(idx,1)-yAmmonia(j));
               min_diff = find(diff==min(diff));
               annotatedIons = [annotatedIons;idx(min_diff)];
           else
               annotatedIons = [annotatedIons;idx];
           end
        else
           continue 
        end
    end
    
    bWater = bseries.waterLoss;
    for j = 1:length(bWater)
        maxDev = ppmDeviation(bWater(j),tolerance);
        idx = find(scanData(:,1) > bWater(j)-maxDev & scanData(:,1) < bWater(j)+maxDev);
        if ~isempty(idx)
           if numel(idx) > 1
               diff = abs(scanData(idx,1)-bWater(j));
               min_diff = find(diff==min(diff));
               annotatedIons = [annotatedIons;idx(min_diff)];
           else
               annotatedIons = [annotatedIons;idx];
           end
        else
           continue 
        end
    end
    
    bAmmonia = yseries.ammoniaLoss;
    for j = 1:length(bAmmonia)
        maxDev = ppmDeviation(bAmmonia(j),tolerance);
        idx = find(scanData(:,1) > bAmmonia(j)-maxDev & scanData(:,1) < bAmmonia(j)+maxDev);
        if ~isempty(idx)
           if numel(idx) > 1
               diff = abs(scanData(idx,1)-bAmmonia(j));
               min_diff = find(diff==min(diff));
               annotatedIons = [annotatedIons;idx(min_diff)];
           else
               annotatedIons = [annotatedIons;idx];
           end
        else
           continue 
        end
    end
    
    filteredSpectrum = scanData(annotatedIons,:);
end