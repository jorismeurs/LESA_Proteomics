classdef quantify
    
    properties
    end
    
    methods
        
        function obj = MS1Quantification(obj)
                        
            reportData = readReport(obj);
            obj.output.reportData = reportData;
            peptideScans = obj.output.reportData(2:end,6);
            peptideAnnotations = obj.output.reportData(2:end,17);
            precursorCharge = obj.output.reportData(2:end,14);
            peptideSequence = obj.output.reportData(2:end,4);

            clc
            for j = 1:length(peptideScans)
               fprintf('(%d) %s | %s \n',j,peptideScans{j},peptideSequence{j}); 
            end
            index = input('Select index for scan of interest: ');
            obj.output.peptideSequence = char(peptideSequence(index));
            obj.output.peptideCharge = char(precursorCharge(index));
            obj.output.peptideCharge = str2num(obj.output.peptideCharge(1:end-1));
            
            if isempty(obj.settings.MS1Tolerance)
                obj.settings.MS1Tolerance = input('MS1 tolerance (ppm): ');
            end
            
            obj = getPeptideMZ(obj,index);
        end
        
        function obj = MS2QuantLibrary(obj)
            % To be developed
            
            % Select peptide
            tempSequence = [];
            for j = 1:length(obj.output.libraryID)
                tempSequence = [tempSequence;{obj.output.libraryID{j}.sequence}'];
            end
            tempSequence = unique(tempSequence);
            clc
            for j = 1:length(tempSequence)
               fprintf('(%d) %s \n',j,tempSequence{j}); 
            end
            peptideIndex = input('Select a peptide: ');
            selectedPeptide = tempSequence{peptideIndex};
            
            % Retrieve best correlation
            for j = 1:length(obj.output.libraryID)
               tempSequence = {obj.output.libraryID{j}.sequence}';
               tempCorr = [obj.output.libraryID{j}.R]';
               matchIndex = find(strcmp(tempSequence,selectedPeptide));
               maxCorr = find(tempCorr(matchIndex,1)==max(tempCorr(matchIndex,1)));
               quantIndex = matchIndex(maxCorr);
               quantSpectrum{j} = obj.output.libraryID{j}(quantIndex).sample;
            end            
            
            % Retrieve top N most intense fragments
            obj.output.topN = [];
            for j = 1:length(quantSpectrum)
                tempSpectrum = cell2mat(quantSpectrum(j));
                [y,b] = fragmentSequence(selectedPeptide);    
                obj.output.topN{j} = retrieveTopNFragments(obj,tempSpectrum,b,y);
            end
        end
        
        function obj = MS2QuantSearch(obj)
            
        end
        
        function obj = getFiles(obj)
            obj.data.fileLocation = [];
            [fileName,pathName] = uigetfile('.mzXML','MultiSelect','on');
            if isequal(fileName,0)
               return 
            end
            obj.data.fileLocation = fullfile(pathName,fileName);
            if ~iscell(obj.data.fileLocation)
                error('Minimum file input is <strong>2</strong>');
            end 
        end
        
        function obj = getMS1Peaks(obj)
            obj.output.MS1Data.totIonCurrent = [];
            obj.data.peakList = [];
            for j = 1:length(obj.data.fileLocation)
               msStruct = mzxmlread(obj.data.fileLocation{j},'Level',1);
               totIonCurrent = [msStruct.scan.totIonCurrent]';
               maxTIC = find(totIonCurrent==(max(totIonCurrent)));
               maxScan = msStruct.scan(maxTIC).peaks.mz;
               obj.output.MS1Data.totIonCurrent(j,1) = msStruct.scan(maxTIC).totIonCurrent;
               tempMZ = maxScan(1:2:end);
               tempInt = maxScan(2:2:end);
               [tempMZ,mzIDX] = unique(tempMZ);
               tempInt = tempInt(mzIDX,1);
               obj.data.peakList{j,1} = mspeaks(tempMZ,tempInt,...
                   'HeightFilter',obj.settings.peakThreshold,...
                   'Denoising',false);
            end 
        end
        
        function obj = getMS1Intensity(obj)
            maxDeviation = ppmDeviation(obj.output.peptideMZ,obj.settings.MS1Tolerance);
            obj.output.MS1Data.peptideIntensity = nan(length(obj.data.fileLocation),1);
            for j = 1:length(obj.data.peakList)
                tempList = cell2mat(obj.data.peakList(j,1));
                matchIndex = find(tempList(:,1) > obj.output.peptideMZ-maxDeviation & ...
                    tempList(:,1) < obj.output.peptideMZ+maxDeviation);
                if ~isempty(matchIndex)
                    if numel(matchIndex) > 1
                        diff = abs(obj.output.peptideMZ-tempList(matchIndex,1));
                        minimumDifference = find(diff==min(diff));
                        obj.output.MS1Data.peptideIntensity(j,1) = tempList(matchIndex(minimumDifference),2);
                    else
                        obj.output.MS1Data.peptideIntensity(j,1) = tempList(matchIndex,2);
                    end
                else
                    continue
                end
            end
        end
        
        function obj = plotResults(obj)
            % Define groups
            
            % Box plot
        end
    end
    
end

