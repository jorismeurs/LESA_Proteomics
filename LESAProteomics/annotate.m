classdef annotate
    
    properties
    end
    
    methods
        function obj = defineParameters(obj)
           clc 
           obj.settings.MS1Tolerance = input('MS1 tolerance (ppm): '); 
           obj.settings.MS2Tolerance = input('MS2 tolerance (Da): ');
        end
        
        function obj = importMGF(obj)
           % Read report data
           reportData = readReport(obj);
            
           % Select .MGF from list
           spectrumFiles = unique(reportData(2:end,5));
           [idx,v] = listdlg('PromptString','Select a spectrum file:',...
               'SelectionMode','single',...
               'ListString',spectrumFiles);
           if v == 0
              return 
           end
           selectedFile = spectrumFiles{idx};
           mgfLocation = fullfile(obj.folder.identification,'data',selectedFile);
           
           % Read .MGF
           obj.output.MGFStruct = readMGF(mgfLocation); 
           obj.output.reportData = reportData;
        end
        
        function obj = annotateMS2(obj)
           
           peptideScans = obj.output.reportData(2:end,6);
           peptideAnnotations = obj.output.reportData(2:end,17);
           precursorCharge = obj.output.reportData(2:end,14);
           peptideSequence = obj.output.reportData(2:end,4);

           clc
           for j = 1:length(peptideScans)
              fprintf('(%d) %s | %s \n',j,peptideScans{j},peptideSequence{j}); 
           end
           scanIndex = input('Select index for scan of interest: ');
           
           
           % Annotate fragment ions
           [yseries,bseries] = fragmentSequence(char(peptideSequence(scanIndex)));
           obj.output.yIons = yseries;
           obj.output.bIons = bseries;
           
           % Make graph
           scanName = peptideScans{scanIndex,1};
           MGFScans = {obj.output.MGFStruct.scan.scanName}';
           MGFIndex = find(strcmp(scanName,MGFScans));
           scanData = obj.output.MGFStruct.scan(MGFIndex).scanData;
           obj = plotMS2Data(obj,scanData)
           
        end
    end
    
end

