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

% Remove PSMs without identification charge
r = [];
r = find(isempty(reportData(:,15)));
if ~isempty(r)
    reportData(r,:) = [];
end

% Remove decoy PSMs
r = [];
r = find(contains(reportData(:,2),'_REVERSED'));
if ~isempty(r)
    reportData(r,:) = [];
end

% Get files and spectrum titles
sequenceList = reportData(2:end,3);
proteinList = reportData(2:end,2);
mgfList = unique(reportData(2:end,10));
scanList = reportData(2:end,11);
precursorCharge = reportData(2:end,15);

% Collect peptide spectra
cd(obj.folder.identification);
count = 0;
tempLibrary  = [];
for j = 1:length(mgfList)
    MGFStruct = readMGF(mgfList{j});
    titleList = {MGFStruct.scan.scanName}';
    for n = 1:length(scanList)
        scanIndex = find(strcmp(titleList,scanList{n}));
        if ~isempty(scanIndex)
           count = count+1; 
           tempLibrary{count,1}  = sequenceList{n};
           tempLibrary{count,2}  = proteinList{n};
           tempLibrary{count,3}  = precursorCharge{n};
           tempLibrary{count,4}  = MGFStruct.scan(scanIndex).scanData;
        end
    end
end

% Retrieve most intense spectrum for unique PSMs
sequence = tempLibrary(:,1);
proteins = tempLibrary(:,2);
charge = tempLibrary(:,3);
scanData = tempLibrary(:,4);

C = unique(sequence);
library = [];
count = 0;
for j = 1:length(C)
    sequenceRow = find(strcmp(sequence,C{j,1}));
    matchCharge = charge(sequenceRow,1);
    charges = unique(matchCharge);
    for n = 1:numel(charges)
        count = count+1;
        chargeIndex = find(strcmp(matchCharge,charges{n}));
        tempTIC = [];
        for k = 1:numel(chargeIndex)
            tempSpectrum = cell2mat(scanData(sequenceRow(chargeIndex,1)));
            tempTIC = sum(tempSpectrum(:,1));
        end
        maxTIC = find(tempTIC==max(tempTIC));
        library(count).sequence = C{j,1};
        library(count).protein = proteins{sequenceRow(1),1};
        library(count).z = charges{n};
        library(count).spectrum = cell2mat(scanData(sequenceRow(chargeIndex(maxTIC),1)));
    end
end

% Copy search files to library folder
cd([obj.folder.identification '\data']);
try
    delete *.cui
    delete *.mzXML
catch
    
end
cd([obj.folder.library '\sample_files']);
delete *.*
movefile('*.mgf',[obj.folder.library '\sample_files']);

% Filter out interfering ions

% Store library as .mat 
obj.output.library = library;
cd([obj.folder.identification '\library']);
save('library.mat','library');
cd(obj.folder.mainFolder);

end

% function filteredSpectrum = retrieveAnnotations(sequence,scanData)
%     [yseries,bseries] = fragmentSequence(sequence);
%     annotatedIons = [];
%     tolerance = obj.settings.MS2Tolerance;
%     
%     yions = yseries.mz;
%     for j = 1:length(yions)
%         idx = find(scanData(:,1) > yions(j)-tolerance & scanData(:,1) < yions(j)+tolerance);
%         if ~isempty(idx)
%            if numel(idx) > 1
%                diff = abs(scanData(idx,1)-yions(j));
%                min_diff = find(diff==min(diff));
%                annotatedIons = [annotatedIons;idx(min_diff)];
%            else
%                annotatedIons = [annotatedIons;idx];
%            end
%         else
%            continue 
%         end
%     end
%     
%     bions = bseries.mz;
%     for j = 1:length(bions)
%         idx = find(scanData(:,1) > bions(j)-tolerance & scanData(:,1) < bions(j)+tolerance);
%         if ~isempty(idx)
%            if numel(idx) > 1
%                diff = abs(scanData(idx,1)-bions(j));
%                min_diff = find(diff==min(diff));
%                annotatedIons = [annotatedIons;idx(min_diff)];
%            else
%                annotatedIons = [annotatedIons;idx];
%            end
%         else
%            continue 
%         end
%     end
%     
%    
%     
%     filteredSpectrum = scanData(annotatedIons,:);
% end