classdef LESAProteomics < visualise & identify
    
    properties (Constant = true)
        version = '0.0.6';
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
            source = fileparts(which('LESA_Proteomics'));           
            addpath(source);
            obj.mainFolder = fullfile(source, 'LESAProteomics');
            addpath(obj.mainFolder);
        end
    end
    
end

