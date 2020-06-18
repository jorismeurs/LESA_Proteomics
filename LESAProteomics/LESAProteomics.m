classdef LESAProteomics < visualise & identify
    
    properties (Constant = true)
        version = '0.1.1';
        developer = 'Joris Meurs, MSc';
    end
    
    properties
        folder
        data
        options
        output
        settings
    end
    
    methods
        function obj = LESAProteomics()
            source = cd;           
            addpath(source);
            obj.folder.mainFolder = fullfile(source);
            obj.folder.identification = fullfile(obj.folder.mainFolder,'identification');
            addpath(obj.folder.identification);
            
            obj.output.MS1Data = [];
            obj.output.MS2Data = [];
            obj.settings.PeptideShakerVersion = '1.16.45';
            obj.settings.SearchGUIVersion = '3.3.20';
        end
    end
    
end

