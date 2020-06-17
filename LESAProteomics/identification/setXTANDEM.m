function settings = setXTANDEM()
    settings.use = '1';
    settings.dynamicRange = input('Dynamic range: ','s');
    settings.minFragmentMZ = input('Minimum fragment m/z: ','s');
    settings.pyro = input('Pyrolidone (0/1)','s');
    settings.acetyl = input('Acetyl (0/1): ','s');
    settings.refinement = input('Refinement (0/1): ','s');
    settings.minPrecMZ = input('Minimum precursor m/z: ','s');
end

