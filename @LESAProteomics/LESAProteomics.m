classdef LESAProteomics 
    
    properties
        options
        output
        parameters
    end
    
    methods
        function obj = LESAProteomics()
            source = fileparts(which('LESAProteomics'));
            source = regexp(source, '.+(?=[@])', 'match'); 
            addpath(source{1});
            
            obj.parameters.
        end
    end
    
end

