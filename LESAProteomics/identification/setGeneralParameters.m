function settings = setGeneralParameters()
    settings.MS1Tolerance = input('MS1 tolerance (ppm): ','s');
    settings.MS2Tolerance = input('MS2 tolerance (Da): ','s');
    settings.miscleavages = input('Miscleavages: ','s');
    settings.minimumCharge = input('Min. charge state: ','s');
    settings.maximumCharge = input('Max. charge state: ','s');
    settings.minPeptideLength = input('Minimum peptide length: ','s');
    settings.maxPeptideLength = input('Maxmimum peptide length: ','s');
    settings.annotationLevel = input('Annotation level (0-100): ','s');
    settings.geneAnnotation = input('Gene annotation (0/1): ','s');
end

