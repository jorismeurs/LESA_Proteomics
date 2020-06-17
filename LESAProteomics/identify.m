classdef identify
    %IDENTIFY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = createParameterFile(obj)
            % Generating a .par file to be used in SearchGUI
            % Files can also be generated in the SearchGUI interface and copied to
            % .../identification/SearchGUI-VERSION

            % Initiate structure
            settings = struct('parameterFileName',[],...
                'xtandem',[],...
                'OMSSA',[],...
                'MSGF',[],...
                'general',[]);

            % Input file name
            clc
            cd([obj.mainFolder '/identification/SearchGUI-3.3.20'])
            settings.parameterFileName = input('Provide file name:   ','s');

            % Enter general parameters
            settings.general = setGeneralParameters();
            clc

            algorithms = {
                'X!Tandem'
                'OMSSA'
                'MSGF+'
                };
            for j = 1:length(algorithms)
                fprintf('(%d) %s \n',j,algorithms{j});
            end
            selectSearchAlgorithms = input('Select search algorithms:   ','s');
            settings.xtandem.use = '0';
            settings.OMSSA.use = '0';
            settings.MSGF.use = '0';

            % Set search algorithm parameters
            if contains(selectSearchAlgorithms,'1')
                settings.xtandem = setXTANDEM();
            elseif contains(selectSearchAlgorithms,'2')
                settings.OMSSA = setOMSSA();
            elseif contains(selectSearchAlgorithms,'3')
                settings.MSGF = setMSGF();
            end
            obj.settings = settings;

            % Generate .par file
            obj = generateFile(settings,obj);

        end
        
        function obj = runID(obj)
            
        end
    end
    
end

