function obj = getPeptideMZ(obj,index)

allMZ = obj.output.reportData(2:end,10);
obj.output.peptideMZ = str2double(allMZ{index});

end

