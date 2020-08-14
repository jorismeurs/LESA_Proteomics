function mz = getMZ(sequence,charge)

z = str2double(charge(1));
aminoAcidMassList = {
    'A' 71.037114
    'R' 156.101111
    'N' 114.042927
    'D' 115.026943
    'C' 103.009185
    'E' 129.042593
    'Q' 128.058578
    'G' 57.021464
    'H' 137.058912
    'I' 113.084064
    'L' 113.084064
    'K' 128.094963
    'M' 131.040485
    'F' 147.068414
    'P' 97.052764
    'S' 87.032028
    'T' 101.047679
    'U' 150.95363
    'W' 186.079313
    'Y' 163.06332	
    'V' 99.068414
    };
mz = [];
for p = 1:numel(sequence)
   idx = find(ismember(aminoAcidMassList(:,1),sequence(p))==true);
   mz = [mz;aminoAcidMassList{idx,2}];
end
mz = sum(mz);

mz = mz+18.010565;

mz = (mz+(z*1.007825)-(z*0.00054858))/z;
disp(mz)


end

