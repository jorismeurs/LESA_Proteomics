classdef LESAProteomics < visualise & identify & annotate
    
    properties (Constant = true)
        version = '0.2.0';
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
            obj.folder.annotation = fullfile(obj.folder.mainFolder,'annotation');
            addpath(obj.folder.identification);
            addpath(obj.folder.annotation);
            
            obj.output.MS1Data = [];
            obj.output.MS2Data = [];
            obj.settings.PeptideShakerVersion = '1.16.45';
            obj.settings.SearchGUIVersion = '3.3.20';
            obj.settings.reportNumber = '13';
        end
    end
    
end

