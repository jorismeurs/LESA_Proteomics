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
            
           % Select PSM from list
           sequenceList = reportData(2:end,3);
           proteinList = reportData(2:end,2);
           mgfList = reportData(2:end,10);
           scanList = reportData(2:end,11);
           clc
           for j = 1:length(proteinList)
               fprintf('(%d) %s | %s \n',j,proteinList{j},sequenceList{j});
           end
           PSMindex = input('Select PSM for annotation:     ');
           mgfLocation = fullfile(obj.folder.identification,'data',mgfList{PSMindex});
           
           % Read .MGF
           obj.output.MGFStruct = readMGF(mgfLocation); 
           obj.output.reportData = reportData;
           obj.output.scanIndex = PSMindex;
        end
        
        function obj = annotateMS2(obj)
           
           peptideScans = obj.output.reportData(2:end,11);
           precursorCharge = obj.output.reportData(2:end,15);
           peptideSequence = obj.output.reportData(2:end,3);
          
           % Annotate fragment ions
           obj.output.peptideSequence = char(peptideSequence(obj.output.scanIndex));
           obj.output.peptideCharge = char(precursorCharge(obj.output.scanIndex));
           [yseries,bseries] = fragmentSequence(char(peptideSequence(obj.output.scanIndex)));
           obj.output.yIons = yseries;
           obj.output.bIons = bseries;
           
           % Make graph
           scanName = peptideScans{obj.output.scanIndex,1};
           MGFScans = {obj.output.MGFStruct.scan.scanName}';
           MGFIndex = find(strcmp(scanName,MGFScans));
           scanData = obj.output.MGFStruct.scan(MGFIndex).scanData;
           obj = plotMS2Data(obj,scanData);
           cd(obj.folder.export);
           saveas(gcf,'');
           cd(obj.folder.mainFolder);
        end
    end
    
end

