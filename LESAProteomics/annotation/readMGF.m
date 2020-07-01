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

% Get scan data and charge (if present)
fprintf('Obtaining MS/MS scans...\n');
beginRow = find(contains(fileData,'BEGIN'));
endRow = find(contains(fileData,'END'));
if length(beginRow) ~= length(endRow)
   error('Unequal rows'); 
end

for n = 1:length(beginRow)
   tempScanData = []; tempMS = []; tempINT = []; tempScan = [];
   tempScanData = fileData(beginRow(n):endRow(n),1);
   chargeRow = find(contains(tempScanData,'CHARGE='));
   if ~isempty(chargeRow)
      tempCharge = fileData{chargeRow,1};
      chargeIDX = find(tempCharge=='='); 
      mgfStruct.scan(n).z = double(str2num(tempCharge(chargeIDX+1:end-1)));
      scanRowStart = chargeRow+1;
      tempMS = str2double(tempScanData(scanRowStart:end-1,1));
      tempMZ = tempMS(1:2:end);
      tempINT = tempMS(2:2:end);
      tempScan = [tempMZ,tempINT];
   else
      mgfStruct.scan(n).z = NaN; 
      scanRowStart = find(contains(tempScanData,'PEPMASS='));
      scanRowStart = scanRowStart+2;
      try
          tempMS = str2double(tempScanData(scanRowStart:end-1,1));
          tempMZ = tempMS(1:2:end);
          tempINT = tempMS(2:2:end);
          tempScan = [tempMZ,tempINT];
      catch
          tempMS = str2double(tempScanData(scanRowStart-1:end-1,1));
          tempMZ = tempMS(1:2:end);
          tempINT = tempMS(2:2:end);
          tempScan = [tempMZ,tempINT];
      end
   end
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



