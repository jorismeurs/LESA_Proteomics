classdef LESAProteomics 
    
    properties (Constant = true)
        version = '0.0.4';
        developer = 'Joris Meurs, MSc';
    end
    
    properties
        mainFolder
        options
        output
        settings
    end
    
    methods
        function obj = LESAProteomics()
            source = fileparts(which('LESAProteomics'));
            source = regexp(source, '.+(?=[@])', 'match');            
            addpath(source{1});
            obj.mainFolder = source{1};
           % addpath(genpath([source{1}, 'annotation']));
           addpath(genpath([source{1}, 'identification']));
           % addpath(genpath([source{1}, 'visualisation']));
           % addpath(genpath([source{1}, 'quantification']));
           % addpath(genpath([source{1}, 'data']));
        end
    end
    
end

