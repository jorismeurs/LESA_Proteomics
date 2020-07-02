function obj = generateFile(obj)
    clc
    searchFolder = [obj.folder.mainFolder '\identification\SearchGUI-' obj.settings.SearchGUIVersion];
    cd(searchFolder);
    jar_struct = dir('*.jar'); 
    searchGUI_jar = jar_struct.name;
    fasta_struct = dir('*.fasta');
    fastaFile = fasta_struct.name;
    generateDecoy = questdlg('Generate decoy sequences?','Decoys',...
        'Yes','No','Yes');
    switch generateDecoy
        case 'Yes'
           system(['cd ' obj.folder.searchFolder]);
           system(['java -cp ' searchGUI_jar ' eu.isas.searchgui.cmd.FastaCLI',...
               ' -in ' fastaFile,...
               ' -decoy']);
           extLoc = find(fastaFile=='.');        
           fastaFile = [fastaFile(1:extLoc-1) '_concatenated_target_decoy'];
        case 'No'
           % Nothing 
        otherwise
           % Handle as 'no' 
    end
    
    
    system(['cd ' obj.folder.searchFolder]);
    system(['java -cp ' searchGUI_jar ' eu.isas.searchgui.cmd.IdentificationParametersCLI',...
        ' -out ' obj.parameters.parameterFileName,...
        ' -db ' fastaFile,...
        ' -prec_tol ' obj.parameters.MS1Tolerance,...
        ' -frag_tol ' obj.parameters.MS2Tolerance,...
        ' -mc ' obj.parameters.miscleavages,...
        ' -min_charge ' obj.parameters.minimumCharge,...
        ' -max_charge ' obj.parameters.maximumCharge,...
        ' -annotation_level ' obj.parameters.annotationLevel,...
        ' -import_peptide_length_min ' obj.parameters.minPeptideLength,...
        ' -import_peptide_length_max ' obj.parameters.maxPeptideLength,...
        ' -useGeneMapping ' obj.parameters.geneAnnotation,...
        ' -db_pi ' fastaFile,...
        ' -xtandem_dynamic_range ' obj.parameters.xtandem.dynamicRange,...
        ' -xtandem_min_frag_mz ' obj.parameters.xtandem.minFragmentMZ,...
        ' -xtandem_quick_acetyl ' obj.parameters.xtandem.acetyl,...
        ' -xtandem_quick_pyro ' obj.parameters.xtandem.pyro,...
        ' -xtandem_refine ' obj.parameters.xtandem.refinement,...
        ' -xtandem_min_prec_mass ' obj.parameters.xtandem.minPrecMZ,...
        ' -msgf_min_pep_length ' obj.parameters.minPeptideLength,...
        ' -msgf_max_pep_length ' obj.parameters.maxPeptideLength,...
        ' -omssa_min_pep_length ' obj.parameters.minPeptideLength,...
        ' -omssa_max_pep_length ' obj.parameters.maxPeptideLength]);
    cd(obj.folder.mainFolder);
end
