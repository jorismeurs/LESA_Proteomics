function noFile = raw2mgf(obj)

[FileName,PathName] = uigetfile({'*.raw','Thermo RAW (.raw)';'*.mgf','MASCOT Generic Format (.mgf)'},...
      'MultiSelect','on');
if isequal(FileName,0)
   noFile = true; 
   return
end

fileLoc = fullfile(PathName,FileName);
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
    for j = 1:length(fileLoc)
        fileInfo = dir(fileLoc{j});
        fileSize = fileInfo.bytes/1024;
        if fileSize > 5000
           idx = [idx;j]; 
        end
    end
else
    idx = 1; 
end

if isequal(fileExt,'raw')
    for j = 1:length(idx)
        system('cd C:\ProteoWizard\');
        system(['msconvert ' fileLoc{idx(j)} ' --mgf --32 --filter "peakPicking true" -o ' [obj.folder.identification '/data']]);
    end
    fileIDX = idx;
    disp('File conversion finished...');
end

% Move .mgf file to data folder
if isequal(fileExt,'mgf')
    copyfile(fileLoc,[obj.folder.identification '\data'])
end

noFile = false;

