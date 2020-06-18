function validateID(obj)

% Add subfolders
addpath([obj.folder.identification '\SearchResults']);
addpath([obj.folder.identification '\PeptideShaker-' obj.settings.PeptideShakerVersion]);
addpath([obj.folder.identification '\SearchGUI-' obj.settings.SearchGUIVersion]);
addpath([obj.folder.identification '\data']);

% Stop any running Java instances
try
   system('taskkill /F /IM java.exe'); 
catch
    
end

% Delete files from existing folders
cd([obj.folder.identification '\data']);
try
    delete *.*
catch 
    warning('No files in data folder');
end
cd([obj.folder.identification '\SearchResults']);
try
    delete *.*
    rmdir('.PeptideShaker_unzip_temp','s');  
catch
    warning('No files in SearchResults folder'); 
end
cd(obj.folder.identification);
try
   delete *.txt
catch
   warning('No .txt files present'); 
end
try
   delete *.cpsx
catch
   warning('No .cpsx files present'); 
end
try
   delete *.html 
catch
   warning('No .html files present'); 
end

cd(obj.folder.mainFolder);

end

