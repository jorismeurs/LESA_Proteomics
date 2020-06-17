classdef LESAProteomics < visualise & identify
    
    properties (Constant = true)
        version = '0.0.7';
        developer = 'Joris Meurs, MSc';
    end
    
    properties
        mainFolder
        data
        options
        output
        settings
    end
    
    methods
        function obj = LESAProteomics()
            source = cd;           
            addpath(source);
            obj.mainFolder = fullfile(source,'LESAProteomics');
            addpath(obj.mainFolder);
            
            obj.output.MS1Data = [];
            obj.output.MS2Data = [];
        end
    end
    
end

