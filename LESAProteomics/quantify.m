classdef quantify
    
    properties
    end
    
    methods
        
        function obj = MS1Quantification(obj)
                        
            reportData = readReport(obj);
            obj.output.reportData = reportData;
            peptideScans = obj.output.reportData(2:end,11);
            precursorCharge = obj.output.reportData(2:end,15);
            peptideSequence = obj.output.reportData(2:end,3);
            mzList = obj.output.reportData(2:end,14);

            clc
            for j = 1:length(peptideScans)
               fprintf('(%d) %s (m/z %.4f) \n',j,peptideSequence{j},str2num(mzList{j})); 
            end
            index = input('Select peptide of interest: ');
            obj.output.peptideSequence = char(peptideSequence(index));
            obj.output.peptideCharge = char(precursorCharge(index));
            obj.output.peptideCharge = str2num(obj.output.peptideCharge(1:end-1));
            
            if isempty(obj.settings.MS1Tolerance)
                obj.settings.MS1Tolerance = input('MS1 tolerance (ppm): ');
            end
            obj.output.peptideMZ = str2double(mzList{index});

        end
        
        function obj = MS2QuantLibrary(obj)           
            % Select peptide
            tempSequence = [];
            for j = 1:length(obj.output.libraryID)
                tempSequence = [tempSequence;cellstr(char(obj.output.libraryID{j}.sequence))];
            end
            tempSequence = unique(tempSequence);
            idx = find(cellfun(@isempty,tempSequence));
            if ~isempty(idx)
                tempSequence(idx,:) = [];
            end
            
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
               if ~isempty(matchIndex)
                    maxCorr = find(tempCorr(matchIndex,1)==max(tempCorr(matchIndex,1)));
                    quantIndex = matchIndex(maxCorr);
                    quantSpectrum{j} = obj.output.libraryID{j}(quantIndex).sample;
               else
                    quantSpectrum{j} = [];
               end
            end            
            
            % Retrieve top N most intense fragments
            obj.output.topN = [];
            for j = 1:length(quantSpectrum)
                tempSpectrum = cell2mat(quantSpectrum(j));
                if ~isempty(tempSpectrum)
                    [y,b] = fragmentSequence(selectedPeptide);    
                    obj.output.topN{j} = retrieveTopNFragments(obj,tempSpectrum,b,y);
                else
                   obj.output.topN{j} = []; 
                end
            end
        end
        
        function obj = MS2QuantSearch(obj)
            % Select peptide
            reportData = readReport(obj);
            obj.output.reportData = reportData;
            decoys = find(contains(reportData(:,2),'_REVERSED'));
            reportData(decoys,:) = [];
            PTMs = find(~cellfun(@isempty,reportData(:,8)));
            reportData(PTMs,:) = [];
            peptideScans = reportData(2:end,11);
            precursorCharge = reportData(2:end,15);
            peptideSequence = reportData(2:end,3);
            mzList = reportData(2:end,14);

            clc
            for j = 1:length(peptideScans)
               fprintf('(%d) %s (m/z %.4f) \n',j,peptideSequence{j},str2num(mzList{j})); 
            end
            index = input('Select peptide of interest: ');
            peptide = peptideSequence{index};
            charge = precursorCharge{index};
            idx = find(strcmp(peptideSequence,peptide) & strcmp(precursorCharge,charge));
            fileData = [reportData(idx+1,10),reportData(idx+1,11)];
            fileData = sortrows(fileData,1);
            obj.output.MS2Data.files = fileData;
            
            % Retrieve MS2 spectrum per file
            obj.output.MS2Data.quantSpectrum = [];
            for j = 1:length(fileData)
                cd([obj.folder.identification '\data']);
                MGFStruct = readMGF(fileData{j,1});
                tempScanNames = {MGFStruct.scan.scanName}';
                scanIndex = find(strcmp(tempScanNames,fileData{j,2}));
                obj.output.MS2Data.quantSpectrum{j,1} = MGFStruct.scan(scanIndex).scanData;
                cd(obj.folder.mainFolder);
            end
            
            % Find matching fragment ions and retrieve top N intensity
            obj.output.MS2Data.topN = [];
            for j = 1:length(obj.output.MS2Data.quantSpectrum)
                tempSpectrum = cell2mat(obj.output.MS2Data.quantSpectrum(j));
                obj.output.MS2Data.totIonCurrent(j) = sum(tempSpectrum(:,2));
                [y,b] = fragmentSequence(peptide);    
                obj.output.MS2Data.topN{j} = retrieveTopNFragments(obj,tempSpectrum,b,y);
            end
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
            if isempty(obj.settings.minMZ)
                obj.settings.minMZ = input('Minimum m/z scan range: ');
            end
            if isempty(obj.settings.maxMZ)
                obj.settings.maxMZ = input('Maximum m/z scan range: ');
            end
            binSize = -8e-8;
            %----------------From Spectral Analysis -----------------
            CMZ = 1/sqrt(obj.settings.minMZ):binSize:1/sqrt(obj.settings.maxMZ)+binSize;
            CMZ = ones(size(CMZ))./(CMZ.^2);
            %--------------------------------------------------------
            
            
            %obj.output.MS1Data.totIonCurrent = [];
            obj.data.peakList = [];
            for j = 1:length(obj.data.fileLocation)
               try
                   msData = mzxml2peaks(mzxmlread(obj.data.fileLocation{j},'Level',1));
                   filteredData = [];
                   for n = 1:length(msData)
                      tempData = cell2mat(msData(n,1));
                      [mz,idx] = unique(tempData(:,1));
                      int = tempData(idx,2);
                      filteredData{n,1} = [mz,int];
                   end
                   interpolatedData = [];
                   fun  =@(x) interp1(x(:,1),x(:,2),CMZ,'linear');
                   interpolatedData = cellfun(fun,filteredData,'UniformOutput',false);
                   interpolatedData = cell2mat(interpolatedData);
                   averageSpectrum = mean(interpolatedData,1);

                   %obj.output.MS1Data.totIonCurrent(j,1) = sum(averageSpectrum);
                   tempMZ = CMZ';
                   tempInt = averageSpectrum';
                   [tempMZ,mzIDX] = unique(tempMZ);
                   tempInt = tempInt(mzIDX,1);
                   obj.data.peakList{j,1} = mspeaks(tempMZ,tempInt,...
                       'HeightFilter',obj.settings.peakThreshold,...
                       'Denoising',false);
               catch
                   obj.data.peakList{j,1} = [];
               end
            end 
        end
        
        function obj = getMS1Intensity(obj)
            maxDeviation = ppmDeviation(obj.output.peptideMZ,obj.settings.MS1Tolerance);
            obj.output.MS1Data.peptideIntensity = nan(length(obj.data.fileLocation),1);
            for j = 1:length(obj.data.peakList)
                tempList = cell2mat(obj.data.peakList(j,1));
                if ~isempty(tempList)
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
                else
                    continue
                end
            end
        end
        
        function obj = libraryPeptideMS1(obj)
            obj.output.MS1Data = [];
            cd(obj.folder.library);
            try
                libraryData = load('library.mat');
            catch 
                error('Library file missing. Run <strong>generateLibraryFile<\strong>');
            end
            cd(obj.folder.mainFolder);
            
            % Select proteins of interest
            allProteinIDs = {libraryData.library.protein}';
            proteinList = unique(allProteinIDs);
            peptideList = {libraryData.library.sequence}';
            chargeList = {libraryData.library.z}';
            selectedProtein = []; k = 0;
            obj.output.selectedLibraryPeptides = [];
            for j = 1:length(proteinList)
               YN = input(sprintf('Include %s ? (Y/N) \n',proteinList{j}),'s');
               if isequal(YN,'Y')                    
                   peptideIDX = find(strcmp(allProteinIDs,proteinList{j,1}));
                   selectedProtein{j,1} = proteinList{j,1};
                   peptides = [];
                   for n = 1:length(peptideIDX)
                        k = k+1;
                        obj.output.selectedLibraryPeptides(k).protein = proteinList{j,1};
                        obj.output.selectedLibraryPeptides(k).peptides = peptideList{peptideIDX(n),1};
                        obj.output.selectedLibraryPeptides(k).z = chargeList(peptideIDX(n),1); 
                   end                  
               else
                   continue
               end
            end
            
            % Calculate peptide m/z for all peptides per protein
            for j = 1:length(obj.output.selectedLibraryPeptides)
               peptideSequence = obj.output.selectedLibraryPeptides(j).peptides;
               charge = obj.output.selectedLibraryPeptides(j).z;
               charge = charge{1};
               [mz,~] = getMZ(peptideSequence,charge);
               obj.output.selectedLibraryPeptides(j).mz = mz;
            end
            
            % Iterate through list and retrieve intensities per protein
            proteinFields = {obj.output.selectedLibraryPeptides.protein}';
            
            for z = 1:length(selectedProtein)                
                fieldNo = find(strcmp(proteinFields,selectedProtein{z}));   
                tempIntensity = [];
                collectIntensity = [];
                for j = 1:length(fieldNo)
                    peptideMZ = obj.output.selectedLibraryPeptides(fieldNo(j)).mz;
                    maxDeviation = ppmDeviation(peptideMZ,obj.settings.MS1Tolerance);
                    for n = 1:length(obj.data.peakList)
                        tempList = cell2mat(obj.data.peakList(n,1));
                        if ~isempty(tempList)
                            matchIndex = find(tempList(:,1) > peptideMZ-maxDeviation & ...
                                tempList(:,1) < peptideMZ+maxDeviation);
                            if ~isempty(matchIndex)
                                if numel(matchIndex) > 1
                                    diff = abs(peptideMZ-tempList(matchIndex,1));
                                    minimumDifference = find(diff==min(diff));
                                    tempIntensity(n,1) = tempList(matchIndex(minimumDifference),2);
                                else
                                    tempIntensity(n,1) = tempList(matchIndex,2);
                                end
                            else
                                tempIntensity(n,1) = NaN;
                                continue
                            end
                        else
                           tempIntensity(n,1) = NaN;
                           continue 
                        end
                    end
                    collectIntensity = [collectIntensity,tempIntensity];
                end
                if isequal(obj.settings.quantType,'average')
                    obj.output.MS1Data(z).proteinAbundance = nanmean(collectIntensity,2);
                    obj.output.MS1Data(z).protein = selectedProtein{z};
                elseif isequal(obj.settings.quantType,'sum')
                    obj.output.MS1Data(z).proteinAbundance = nansum(collectIntensity,2);
                    obj.output.MS1Data(z).protein = selectedProtein{z};
                else
                    error('Quantification setting should be either <strong>average</strong> or <strong>sum</strong>');
                end
            end
        end
        
        function obj = setGroups(obj)
            reportData = readReport(obj);
            fileList = unique(reportData(2:end,10));

            obj.settings.group_code = [];
            for j = 1:length(fileList)
                obj.settings.group_code = [obj.settings.group_code;input(sprintf('Group no. for %s:  ',fileList{j}))];                
            end
            obj.settings.groups = unique(obj.settings.group_code);
            
            obj.settings.groupName = [];
            for j = 1:length(obj.settings.groups)
                obj.settings.groupName = [obj.settings.groupName;cellstr(input(sprintf('Group label for %s: ',num2str(obj.settings.groups(j))),'s'))];
            end
        end
        
        function obj = plotResultsMS1(obj)
                       
            data = obj.output.MS1Data.peptideIntensity;
            
            % Jitter plot
            jitterplot(log(data),obj.settings.group_code,...
                'meanWidth',0.06,...
                'stdWidth',0.1)
            set(gca,'XTick',1:1:length(obj.settings.groups));
            set(gca,'XTickLabel',obj.settings.groupName);
            set(gca,'FontName','Calibri','FontSize',14);
            set(gcf,'Color','white');
            ylabel('Log intensity');
            
            cd(obj.folder.export);
            saveas(gcf,['MS1_Quant' obj.settings.imageFormat]);
            cd(obj.folder.mainFolder);
        end
        
        function obj = plotResultsMS2(obj)
                       
            data = [];
            for j = 1:length(obj.output.MS2Data.topN)
                tempData = cell2mat(obj.output.MS2Data.topN(j));
                data = [data;max(tempData)];
            end
            
            % Jitter plot
            jitterplot(log(data),obj.settings.group_code,...
                'meanWidth',0.06,...
                'stdWidth',0.1)
            set(gca,'XTick',1:1:length(obj.settings.groups));
            set(gca,'XTickLabel',obj.settings.groupName);
            set(gca,'FontName','Calibri','FontSize',14);
            set(gcf,'Color','white');
            ylabel('Log intensity');
            
            cd(obj.folder.export);
            saveas(gcf,['MS2_Quant' obj.settings.imageFormat]);
            cd(obj.folder.mainFolder);
        end
        
        function obj = peptideMassFingerprinting(obj)
            
            if isempty(obj.settings.minMZ)
                obj.settings.minMZ = input('Minimum m/z scan range: ');
            end
            if isempty(obj.settings.maxMZ)
                obj.settings.maxMZ = input('Maximum m/z scan range: ');
            end
            binSize = -8e-8;
            %----------------From Spectral Analysis -----------------
            CMZ = 1/sqrt(obj.settings.minMZ):binSize:1/sqrt(obj.settings.maxMZ)+binSize;
            CMZ = ones(size(CMZ))./(CMZ.^2);
            %-------------------------------------------------------- 
           
            % Load exported peptide list generated in mMass
            [fileName,pathName] = uigetfile('.xlsx');
            if isequal(fileName,0)
                disp('No file selected');
                return
            end
            peptideFile = fullfile(pathName,fileName);
            
            % Load mzXML files of all samples
            [fileName,pathName] = uigetfile('.mzXML','MultiSelect','on');
            if isequal(fileName,0)
                disp('No file selected');
                return
            end
            mzXMLFiles = fullfile(pathName,fileName);
            
            % Get peak list for each file based on set parameters
            for j = 1:length(mzXMLFiles)
                scanList = mzxml2peaks(mzxmlread(mzXMLFiles{j},'Level',1));
                
                % Generate average MS1 spectrum               
                interpolatedData = []; averageSpectrum = [];
                for n = 1:length(scanList)
                    tempData = cell2mat(scanList(n,1));
                    [mz,idx] = unique(tempData(:,1));
                    int = tempData(idx,2);
                    filteredData{n,1} = [mz,int];
                end              
                fun  =@(x) interp1(x(:,1),x(:,2),CMZ,'linear');
                interpolatedData = cellfun(fun,filteredData,'UniformOutput',false);
                interpolatedData = cell2mat(interpolatedData);
                averageSpectrum = mean(interpolatedData,1);
                
                peakList{j,1} = mspeaks(CMZ,averageSpectrum,'HeightFilter',obj.settings.peakThreshold);
            end
            
            % Match theoretical m/z values with observed values
            [~,sheets] = xlsfinfo(peptideFile);
            for z = 1:numel(sheets)
                [peptideData,peptideTXT] = xlsread(peptideFile,sheets{z});
                for j = 1:length(peakList)
                   tempPeaks = cell2mat(peakList(j,1));
                   tempPeptideInt = [];
                   for n = 1:size(peptideData,1)
%                        if tempPeaks(n,1) < obj.settings.minMZ || tempPeaks(n,1) > obj.settings.maxMZ
%                            continue
%                        end
                       idx = find(tempPeaks(:,1) > peptideData(n,2)-obj.settings.MS1Tolerance & ...
                           tempPeaks(:,1) < peptideData(n,2)+obj.settings.MS1Tolerance);
                       if ~isempty(idx)
                           if numel(idx) > 1
                               diff = abs(tempPeaks(idx,1)-peptideData(n,2));
                               min_diff = find(diff==min(diff));
                               tempPeptideInt = [tempPeptideInt;tempPeaks(idx(min_diff),2)];
                           else
                              tempPeptideInt = [tempPeptideInt;tempPeaks(idx,2)]; 
                           end
                       else
                          tempPeptideInt = [tempPeptideInt;0]; 
                       end                           
                   end
                   obj.output.PMF(z).protein(j).FileName = fileName{j};
                   obj.output.PMF(z).protein(j).proteinAbundance = sum(tempPeptideInt);
                end
                
            end
        end
    end
    
end

