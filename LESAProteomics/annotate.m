classdef annotate
    
    properties
    end
    
    methods
        function obj = defineParameters(obj)
           clc 
           obj.settings.MS1Tolerance = input('MS1 tolerance (ppm): '); 
           obj.settings.MS2Tolerance = input('MS2 tolerance (Da): ');
        end
        
        function obj = selectMS1(obj)
            clc
            reportData = readReport(obj);
            obj.output.reportData = reportData;
            proteinList = [];  
            proteinList = reportData(2:end,2);
            decoys = find(contains(proteinList,'_REVERSED'));
            decoys = decoys+1;
            reportData(decoys,:) = [];
            mzList = []; proteinList = []; chargeList = [];
            mzList = reportData(2:end,14);
            proteinList = reportData(2:end,2);
            chargeList = reportData(2:end,15);
                        
            % Select file
            fileList = unique(reportData(2:end,10));
            allFiles = reportData(2:end,10);
            for j = 1:length(fileList)
                fprintf('(%d) %s \n',j,fileList{j});
            end
            fileIndex = input('Choose file for annotation:  ');
            mzIndex = find(strcmp(allFiles,fileList{fileIndex}));

            obj.output.mzList = mzList(mzIndex,1);
            obj.output.proteinList = proteinList(mzIndex,1);
            obj.output.chargeList = chargeList(mzIndex,1);
            fileType = questdlg('File type','Select file type',...
                '.RAW','.mzXML','.RAW');
            disp(fileType);
            switch fileType
                case '.RAW'
                    raw2mzxml(obj);
                case '.mzXML'
                    selectMZXML(obj);
                otherwise
                    return
            end
            
            fileName = fileList{fileIndex};
            extensionIndex = find(fileName=='.');
            obj.output.MS1File = fileName(1:extensionIndex-1);
            obj.output.mzXMLFile = [obj.folder.identification '\data\' fileName(1:extensionIndex-1) '.mzXML'];
            if isequal(fileType,'.RAW')
                scanPeakList = mzxml2peaks(mzxmlread(obj.output.mzXMLFile));
                bpIntensity = cellfun(@(x) max(x(:,2)),scanPeakList);
                maxIndex = find(bpIntensity==max(bpIntensity));  
                obj.output.scanData = cell2mat(scanPeakList(maxIndex,1));
            elseif isequal(fileType,'.mzXML')
                scanPeakList = cell2mat(mzxml2peaks(mzxmlread(obj.output.mzXMLFile)));
                obj.output.scanData = mspeaks(scanPeakList(:,1),scanPeakList(:,2));
            end
        end
        
        function obj = annotateMS1(obj)           
            plotMS1Data(obj);
            cd(obj.folder.export);
            saveas(gcf,[obj.output.MS1File obj.settings.imageFormat]);
            cd(obj.folder.mainFolder);
        end
        
        function obj = annotateMS2(obj)
           % Read report data
           reportData = readReport(obj);
           
           % Filter PSM
           scores = reportData(2:end,24);
           numericScores = [];
           for j = 1:length(scores)
              numericScores = [numericScores;str2num(char(scores(j,1)))];
           end
           removeIDX = find(numericScores(:,1)<obj.settings.minPSMScore | isnan(numericScores(:,1)));
           removeIDX = removeIDX+1;
           reportData(removeIDX,:) = [];
            
           % Select PSM from list
           sequenceList = reportData(2:end,3);
           proteinList = reportData(2:end,2);
           mgfList = reportData(2:end,10);
           scanList = reportData(2:end,11);
           mzList = reportData(2:end,14);
           precursorCharge = reportData(2:end,16);
           clc
           for j = 1:length(proteinList)
               fprintf('(%d) %s | %s (m/z %s %s) \n',j,proteinList{j},sequenceList{j},mzList{j},precursorCharge{j});
           end
           PSMindex = input('Select PSM for annotation:     ');
           mgfLocation = fullfile(obj.folder.identification,'data',mgfList{PSMindex});
           
           % Read .MGF
           obj.output.MGFStruct = readMGF(mgfLocation); 
           obj.output.reportData = reportData;
           obj.output.scanIndex = PSMindex;
           
           peptideScans = obj.output.reportData(2:end,11);
           precursorCharge = obj.output.reportData(2:end,16);
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
           saveas(gcf,[obj.output.peptideSequence obj.settings.imageFormat]);
           cd(obj.folder.mainFolder);
        end
        
        function obj = annotateLibrary(obj)
           clc
           fileNo = input(sprintf('Select file (1-%d): ',length(obj.output.libraryID)));
           tempIDs = obj.output.libraryID{fileNo};
           tempProtein = {tempIDs.protein}';
           tempSequence = {tempIDs.sequence}';
           tempCharge = {tempIDs.z}';
           for j = 1:length(tempProtein)
               fprintf('(%d) %s | %s \n',j,tempProtein{j,1},tempSequence{j,1});
           end
           matchIndex = input('Select match for annotation:     ');
           
           % Plot data
           obj = plotLibraryData(obj,tempIDs(matchIndex).reference,tempIDs(matchIndex).sample,tempIDs(matchIndex).originalSpectrum,tempIDs(matchIndex).fragmentIons,tempSequence{matchIndex},tempCharge{matchIndex},tempIDs(matchIndex).R);
           cd(obj.folder.export);
           saveas(gcf,[tempSequence{matchIndex,1} '_library' obj.settings.imageFormat]);
           cd(obj.folder.mainFolder);
        end
    end
    
end

