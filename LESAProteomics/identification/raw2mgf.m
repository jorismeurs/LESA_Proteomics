function noFile = raw2mgf(outputFolder)

[FileName,PathName] = uigetfile({'*.raw','Thermo RAW (.raw)';'*.mgf','MASCOT Generic Format (.mgf)'},...
      'MultiSelect','on');
if isequal(FileName,0)
   noFile = true; 
   return
end

fileLoc = fullfile(PathName,FileName);
% spaceIDX = find(fileLoc==' ');
% if ~isempty(spaceIDX)
%    error('File path contains spaces'); 
% end

if iscell(fileLoc)
    fileExtLoc = find(FileName{1}=='.');
    fileExt = FileName{1}(fileExtLoc+1:end);
else
    fileExtLoc = find(FileName=='.');
    fileExt = FileName(fileExtLoc+1:end);
end

% If file size < 2 MB, skip file
clc
idx = [];
if iscell(fileLoc)
    idx = 1:1:length(fileLoc); 
else
    idx = 1; 
end

if isequal(fileExt,'raw')
    for j = 1:length(idx)
        system('cd C:\ProteoWizard\');
        system(['msconvert "' fileLoc{idx(j)} '" --mgf --32 --filter "peakPicking true" -o ' outputFolder]);
    end
    fileIDX = idx;
    disp('File conversion finished...');
end

% Move .mgf file to data folder
if isequal(fileExt,'mgf')
    for j = 1:length(fileLoc)
        copyfile(fileLoc{j},outputFolder)
    end
end

noFile = false;

