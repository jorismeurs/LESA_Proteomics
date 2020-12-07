function selectMZXML(obj)

[fileName,pathName] = uigetfile('.mzXML',...
    'Select .mzXML file');
if isequal(fileName,0)
    return
end
fileLocation = fullfile(pathName,fileName);
%copyfile(fileLocation,[obj.folder.identification '/data'])


end

