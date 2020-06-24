function obj = results2library(obj)

% Open report file
reportData = readReport(obj);

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
           library(count).spectrum = MGFStruct.scan(scanIndex).scanData;
        end
    end
end
obj.output.library = library;
cd([obj.folder.identification '\library']);
save('library.mat','library');
cd(obj.folder.mainFolder);

end

