classdef identify
    
    properties
    end
    
    methods
        function obj = createParameterFile(obj)
            % Generating a .par file to be used in SearchGUI
            % Files can also be generated in the SearchGUI interface and copied to
            % .../identification/SearchGUI-VERSION
            if ~isempty(obj.parameters) && any(structfun(@isempty,obj.parameters)) 
                executeFileGeneration = questdlg(...
                    'Use current input as parameters?',...
                    'Generate parameter file',...
                    'Yes','No','Yes');
            else
                executeFileGeneration = 'Yes';
            end
            
            switch executeFileGeneration
                case 'No'
                    clc
                    obj.folder.searchFolder = [obj.folder.identification '/SearchGUI-' obj.settings.SearchGUIVersion];
                    cd(obj.folder.searchFolder);
                    obj.parameters.parameterFileName = input('Provide file name:   ','s');

                    obj = setGeneralParameters(obj);
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
                    obj.parameters.xtandem.use = '0';
                    obj.parameters.OMSSA.use = '0';
                    obj.parameters.MSGF.use = '0';

                    if contains(selectSearchAlgorithms,'1')
                        obj.parameters.xtandem = setXTANDEM();
                    elseif contains(selectSearchAlgorithms,'2')
                        obj.parameters.OMSSA = setOMSSA();
                    elseif contains(selectSearchAlgorithms,'3')
                        obj.parameters.MSGF = setMSGF();
                    end

                    generateFile(obj);
                case 'Yes'
                    generateFile(obj);
            end
        end
        
        function obj = runSearch(obj)
            validateID(obj);
            noFile = raw2mgf(obj);
            if noFile == true
                error('No files selected');
            end
            
            searchGUILocation = [obj.folder.identification '\SearchGUI-' obj.settings.SearchGUIVersion];
            paramLocation = searchGUILocation;
            cd(searchGUILocation);
            if isempty(dir('*.par'))
               disp('No parameter file present')
               return 
            else
               parStruct = dir('*.par'); 
               obj.settings.parameterFileName = fullfile(parStruct.folder,parStruct.name);
            end
            
            try
                searchEngines = [' -xtandem ' obj.settings.xtandem.use ' -omssa ' obj.settings.OMSSA.use ' -msgf ' obj.settings.MSGF.use];    
            catch
                options = {
                    'xtandem'
                    'omssa'
                    'msgf'
                    };
                clc
                for j = 1:length(options)
                    fprintf('(%d) %s \n',j,options{j});
                end
                searchInput = input('Select search algorithm: ','s');
                if contains(searchInput,'1')
                    obj.settings.xtandem.use = '1';
                else
                    obj.settings.xtandem.use = '0';
                end
                if contains(searchInput,'2')
                    obj.settings.OMSSA.use = '1';
                else
                    obj.settings.OMSSA.use = '0';
                end
                if contains(searchInput,'3')
                    obj.settings.MSGF.use = '1';
                else
                    obj.settings.MSGF.use = '0';
                end
                searchEngines = [' -xtandem ' obj.settings.xtandem.use ' -omssa ' obj.settings.OMSSA.use ' -msgf ' obj.settings.MSGF.use];
            end
                        
            outputFolder = [obj.folder.identification '\SearchResults'];
            inputFolder = [obj.folder.identification '\data'];

            cd(searchGUILocation);
            jar_struct = dir('*.jar'); 
            searchGUI_jar = jar_struct.name;

            system(['cd ' searchGUILocation],'-echo');
            system(['java -cp ' searchGUI_jar ' eu.isas.searchgui.cmd.SearchCLI -spectrum_files ' inputFolder ' -output_folder ' outputFolder ' -id_params ' obj.settings.parameterFileName searchEngines ' -output_option 1']);
            
            cd(obj.folder.mainFolder);
        end
        
        function obj = PeptideShaker(obj)
            clc

            peptideShakerLocation = [obj.folder.identification '\PeptideShaker-' obj.settings.PeptideShakerVersion];
            spectrumFiles = [obj.folder.identification '\data'];
            identificationFiles = [obj.folder.identification '\SearchResults'];

            try
                cd(peptideShakerLocation);
                jar_struct = dir('*.jar'); 
                peptideShaker_jar = jar_struct.name;
            catch
                error('PeptideShaker .jar file not found'); 
            end

            outputLabel = ['output_' datestr(datetime('now'),'yyyymmdd_hhMMss') '.cpsx'];
            system(['cd ' peptideShakerLocation],'-echo');
            system(['java -Xms1024m -Xmx1024m -cp ' peptideShaker_jar ' eu.isas.peptideshaker.cmd.PeptideShakerCLI -experiment A_ -sample B_ -replicate 1 -identification_files ' identificationFiles ' -spectrum_files ' spectrumFiles ' -out ' fullfile(obj.folder.identification,outputLabel)]);
            obj.output.PeptideShaker = fullfile(obj.folder.identification,outputLabel);
            
            cd(obj.folder.mainFolder);
        end
        
        function obj = generateReport(obj)
            clc

            % Execute ReportCLI
            cd(obj.folder.identification);
            try
                cpsx_struct = dir('*.cpsx');
                cpsxLocation = fullfile(cpsx_struct.folder,cpsx_struct.name);
            catch
                cd(obj.folder.mainFolder);
                error(sprintf('No .cpsx file found in %s',obj.folder.identification)); 
            end
            
            peptideShakerLocation = [obj.folder.identification '\PeptideShaker-' obj.settings.PeptideShakerVersion];
            try
                cd(peptideShakerLocation);
                jar_struct = dir('*.jar'); 
                peptideShaker_jar = jar_struct.name;
            catch
                error('PeptideShaker .jar file not found'); 
            end

            cd(peptideShakerLocation);
            system(['cd ' peptideShakerLocation],'-echo');
            system(['java -Xms1024m -Xmx1024m -cp ' peptideShaker_jar ' eu.isas.peptideshaker.cmd.ReportCLI -in ' cpsxLocation ' -out_reports ' obj.folder.identification ' -reports ' obj.settings.reportNumber]); 
            
            cd(obj.folder.mainFolder);
        end
        
        function obj = generateLibraryFile(obj)
            obj = results2library(obj);
        end
        
        function obj = libraryIdentification(obj)
            raw2mgf(obj);
            obj = getSampleSpectra(obj);
        end
    end
    
end

