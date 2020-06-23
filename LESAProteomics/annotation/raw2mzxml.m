function raw2mzxml(obj)

[fileName,pathName] = uigetfile('.RAW',...
    'Select .RAW file for conversion');
if isequal(fileName,0)
    return
end
fileLocation = fullfile(pathName,fileName);

system('cd C:\ProteoWizard\');
system(['msconvert "' fileLocation '" --mzXML --32 --filter "peakPicking true" "msLevel 1" -o ' [obj.folder.identification '/data']]);


end

