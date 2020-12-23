classdef LESAProteomics < visualise & identify & annotate & quantify
    
    properties (Constant = true)
        version = '0.10.2';
        developer = 'Joris Meurs, MSc';
    end
    
    properties
        folder
        data
        output
        settings
        parameters
    end
    
    methods
        function obj = LESAProteomics()
            source = fileparts(which('LESAProteomics'));           
            addpath(source);
            cd(source);
            obj.folder.mainFolder = fullfile(source);
            obj.folder.identification = fullfile(obj.folder.mainFolder,'identification');
            obj.folder.annotation = fullfile(obj.folder.mainFolder,'annotation');
            obj.folder.quantification = fullfile(obj.folder.mainFolder,'quantification');
            obj.folder.export = fullfile(obj.folder.mainFolder,'images');
            obj.folder.library = fullfile(obj.folder.identification,'library');
            addpath(obj.folder.identification);
            addpath(obj.folder.annotation);
            addpath(obj.folder.quantification);
            addpath(obj.folder.export);
            addpath(obj.folder.library);
            
            obj.output.MS1Data = [];
            obj.output.MS2Data = [];
            obj.output.annotations = [];
            obj.output.file = [];
            obj.output.file.identification = [];
            obj.output.selectedLibraryPeptides = [];
            obj.settings.PeptideShakerVersion = '1.16.45';
            obj.settings.SearchGUIVersion = '3.3.20';
            obj.settings.reportNumber = '11';
            obj.settings.minPSMScore = 95;
            obj.settings.topNFragments = 3;
            obj.settings.XLim = [];
            obj.settings.YLim = [];
            obj.settings.fontSize = 11;
            obj.settings.minMZ = [];
            obj.settings.maxMZ = [];
            obj.settings.minCharge = 2;
            obj.settings.maxCharge = 4;
            
            obj.settings.precursorRemoval = true;
            obj.settings.MS1Tolerance = 10; % ppm
            obj.settings.MS2Tolerance = 0.02; % Da
            obj.settings.peakThreshold = 1e4;
            obj.settings.imageFormat = '.tif';
            obj.settings.neutralLoss = false;
            
            obj.settings.fontSize = 9; % pt
            obj.settings.width = 1; % pt
            
            obj.folder.searchFolder = [obj.folder.identification '/SearchGUI-' obj.settings.SearchGUIVersion];
        end
    end
    
end

