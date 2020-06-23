function obj = setGeneralParameters(obj)
    obj.parameters.MS1Tolerance = input('MS1 tolerance (ppm): ','s');
    obj.parameters.MS2Tolerance = input('MS2 tolerance (Da): ','s');
    obj.parameters.miscleavages = input('Miscleavages: ','s');
    obj.parameters.minimumCharge = input('Min. charge state: ','s');
    obj.parameters.maximumCharge = input('Max. charge state: ','s');
    obj.parameters.minPeptideLength = input('Minimum peptide length: ','s');
    obj.parameters.maxPeptideLength = input('Maxmimum peptide length: ','s');
    obj.parameters.annotationLevel = input('Annotation level (0-100): ','s');
    obj.parameters.geneAnnotation = input('Gene annotation (0/1): ','s');
end

