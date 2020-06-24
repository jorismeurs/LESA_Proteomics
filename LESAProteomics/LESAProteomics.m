classdef LESAProteomics < visualise & identify & annotate & quantify
    
    properties (Constant = true)
        version = '0.5.0';
        developer = 'Joris Meurs, MSc';
    end
    
    properties
        folder
        data
        options
        output
        settings
        parameters
    end
    
    methods
        function obj = LESAProteomics()
            source = cd;           
            addpath(source);
            obj.folder.mainFolder = fullfile(source);
            obj.folder.identification = fullfile(obj.folder.mainFolder,'identification');
            obj.folder.annotation = fullfile(obj.folder.mainFolder,'annotation');
            obj.folder.quantification = fullfile(obj.folder.mainFolder,'quantification');
            obj.folder.export = fullfile(obj.folder.mainFolder,'images');
            addpath(obj.folder.identification);
            addpath(obj.folder.annotation);
            addpath(obj.folder.quantification);
            addpath(obj.folder.export);
            
            obj.output.MS1Data = [];
            obj.output.MS2Data = [];
            obj.settings.PeptideShakerVersion = '1.16.45';
            obj.settings.SearchGUIVersion = '3.3.20';
            obj.settings.reportNumber = '11';
            
            obj.settings.MS1Tolerance = [];
            obj.settings.MS2Tolerance = [];
            obj.settings.peakThreshold = 1e4;
            obj.settings.imageFormat = '.tif';
            
            obj.folder.searchFolder = [obj.folder.identification '/SearchGUI-' obj.settings.SearchGUIVersion];
        end
    end
    
end

