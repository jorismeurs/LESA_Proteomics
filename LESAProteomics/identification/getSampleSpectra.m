function obj = getSampleSpectra(obj)

[fileName,pathName] = uigetfile('.mgf','MultiSelect','on');
if isequal(fileName,0)
    return
end
fileLocation = fullfile(pathName,fileName);

if ~iscell(fileLocation)
   obj.output.sampleMGF = readMGF(fileLocation);
else
   for j = 1:length(fileLocation)
       obj.output.sampleMGF{j} = readMGF(fileLocation{j});
   end
end


end

