function obj = createParameterFile(obj)
% Generating a .par file to be used in SearchGUI
% Files can also be generated in the SearchGUI interface and copied to
% /@LESAProteomics/identification/SearchGUI-VERSION

% Input file name
clc
parameterFileName = input('Provide file name:   ','s');

% Enter general parameters
settings = setGeneralParameters();
clc

algorithms = {
    'X!Tandem'
    'OMSSA'
    'MSGF+'
    };
for j = 1:length(algorithms)
    fprintf('(%d) %s',j,algorithms{j});
end
selectSearchAlgorithms = input('Select search algorithms:   ','s');
settings.xtandem.use = false;

% Set search algorithm parameters
if contains(selectSearchAlgorithms,'1')
    settings = setXTANDEM();
elseif contains(selectSearchAlgorithms,'2')
    settings = setOMSSA();
elseif contains(selectSearchAlgorithms,'3')
    settings = setMSGF();
end

% Generate .par file
generateFile(parameterFileName,settings)

end

function settings = setGeneralParameters()
    settings.MS1Tolerance = input('MS1 tolerance (ppm): ');
    settings.MS2Tolerance = input('MS2 tolerance (Da): ');
    settings.miscleavages = input('Miscleavages: ');
    settings.minimumCharge = input('Min. charge state: ');
    settings.maximumCharge = input('Max. charge state: ');
    settings.minPeptideLength = input('Minimum peptide length: ');
    settings.maxPeptideLength = input('Maxmimum peptide length: ');
end

function settings = setXTANDEM()
    settings.xtandem.use = true;
    settings.xtandem.dynamicRange = input('Dynamic range: ');
    settings.xtandem.minFragmentMZ = input('Minimum fragment m/z: ');
    settings.xtandem.pyro = input('Pyrolidone (0/1)');
    settings.xtandem.acetyl = input('Acetyl (0/1): ');
    settings.xtandem.refinement = input('Refinement (0/1): ');
    settings.xtandem.minPrecMZ = input('Minimum precursor m/z: ');
end

function settings = setOMSSA()
    settings.OMSSA.use = true;
end

function settings = setMSGF()
    settings.MSGF.use = true;
end