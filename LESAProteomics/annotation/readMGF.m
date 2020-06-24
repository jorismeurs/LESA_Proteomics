function mgfStruct = readMGF(file)

% Validate file
ismgf(file);

% Default structure for data storage
defaultStruct = struct('scanName',[],...
    'precursorMass',[],...
    'precursorIntensity',[],...
    'z',[],...
    'scanData',[]);

fprintf('Reading file\n');
% Open and read .MGF file

fileID = fopen(file,'r');   
fileData = [];
fileData = textscan(fileID,'%s');
fileData = fileData{1,1};
fclose(fileID);

% Get total number of scans
scanCount = numel(find(contains(fileData,'BEGIN')));
scan = repmat(defaultStruct,scanCount,1);
fprintf('Number of MS/MS scans: %d\n',scanCount);

% Get scan names
fprintf('Obtaining scan titles...\n');
try
    titleRow = find(contains(fileData,'TITLE='));
    for n = 1:length(titleRow)
       tempTitle = []; titleIDX = [];
       tempTitle = fileData{titleRow(n),1};
       titleIDX = find(tempTitle=='=');
       mgfStruct.scan(n).scanName = tempTitle(titleIDX+1:end);
    end
catch
    warning('Titles not found');
end

% Get precursor mass & intensity
fprintf('Obtaining precursor data...\n');
precursorRow = find(contains(fileData,'PEPMASS='));
try
    for n = 1:length(precursorRow)
       tempPrecursor = []; precursorIDX = []; precursorInt = [];
       tempPrecursor = fileData{precursorRow(n),1};
       precursorIDX = find(tempPrecursor=='=');
       mgfStruct.scan(n).precursorMass = double(str2num(tempPrecursor(precursorIDX+1:end)));
       precursorInt = double(str2num(fileData{precursorRow(n)+1,1}));
       mgfStruct.scan(n).precursorIntensity = precursorInt;
    end
catch 
    warning('Precursor data not found')
end

% Get charge
fprintf('Obtaining precursor charge...\n');
chargeRow = find(contains(fileData,'CHARGE='));
try
    for n = 1:length(chargeRow)
       tempCharge = []; chargeIDX = [];
       tempCharge = fileData{chargeRow(n),1};
       chargeIDX = find(tempCharge=='=');
       mgfStruct.scan(n).z = double(str2num(tempCharge(chargeIDX+1:end-1)));
    end
catch
    warning('Charge not found');
end

% Get MS scan data (m/z and intensities)
fprintf('Obtaining MS/MS scans...\n');
startRow = find(contains(fileData,'CHARGE='));
if isempty(startRow)
    warning('Charge state not found');
    startRow = find(contains(fileData,'PEPMASS='));
    startRow = startRow+2; % Second row after 'PEPMASS' is first m/z value
else
    startRow = startRow+1; % Spectral data starts in row after charge
end
disp(startRow(1));
endRow = find(contains(fileData,'END'));
endRow = endRow-1; % Spectral data ends row before 'END'
if length(startRow) ~= length(endRow)
   error('Unequal rows'); 
end
if isempty(startRow) || isempty(endRow)
   warning('No MS/MS data');
end
for n = 1:length(startRow)
    tempMS = []; tempMZ = []; tempINT = []; tempScan = [];
    tempMS = str2double(fileData(startRow(n):endRow(n),1));
    tempMZ = tempMS(1:2:end);
    tempINT = tempMS(2:2:end);
    tempScan = [tempMZ,tempINT];
    mgfStruct.scan(n).scanData = tempScan;
end

fprintf('Finished!\n');

end

function ismgf(files)   
    if ~iscell(files)
        tf = ~isempty(find(contains(files,'mgf')));
        if tf == true
            fprintf('File type: .mgf\n');
        else
            error('File type is not .mgf');
        end
    else
        for p = 1:length(files)
            tf = ~isempty(find(contains(files{p},'mgf')));
            if tf == true
                fprintf('File %d. Type is .mgf\n',p);
            else
                error(sprintf('File %d. Type is not .mgf',p));
            end
        end
    end  
end



