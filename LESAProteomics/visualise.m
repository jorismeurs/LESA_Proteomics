classdef visualise
   
    properties
        
    end
    
    methods
        function obj = importData(obj)            
            [FileName,PathName] = uigetfile('.mzXML','MultiSelect','on');
            if isequal(FileName,0)
               return 
            end
            filePath = fullfile(PathName,FileName);

            if ~iscell(filePath)
               fileCount = 1; 
            else
               fileCount = length(filePath); 
            end

            for j = 1:fileCount
                if fileCount == 1
                    msStruct = mzxmlread(filePath,'Level',1);
                else
                    msStruct = mzxmlread(filePath{j},'Level',1);
                end
                totalIonCount = [msStruct.scan.totIonCurrent];
                maxTIC = find(totalIonCount==max(totalIonCount));
                data{j} = msStruct.scan(maxTIC).peaks.mz;
            end
            obj.data.MS1Data = data;
        end
        
        function plotMS1(obj)

            hold on
            for j = 1:length(obj.data.MS1Data)
                tempData = cell2mat(obj.data.MS1Data(j));
                mz = tempData(1:2:end);
                int = tempData(2:2:end);
                plot(mz,int);
            end
            xlim([min(mz) max(mz)]);
            hold off
            xlabel('m/z');
            ylabel('Intensity');
            set(gca,'FontName','Calibri');
            set(gcf,'Color','white');
        end
    end
    
end