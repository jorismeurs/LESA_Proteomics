function topN = retrieveTopNFragments(obj,spectrum,b,y)
%
% Most intense b and y ions

YIons = [y.mz]';
BIons = [b.mz]';

intensities = [];
for j = 1:length(YIons)
    matchIon = find(spectrum(:,1) > YIons(j)-obj.settings.MS2Tolerance & ...
        spectrum(:,1) < YIons(j)+obj.settings.MS2Tolerance);
    if ~isempty(matchIon)
       if numel(matchIon) > 1
          diff = abs(YIons(j)-spectrum(matchIon,1));
          min_diff = find(diff==min(diff));
          intensities = [intensities;spectrum(matchIon(min_diff),2)];
       else
          intensities = [intensities;spectrum(matchIon,2)];
       end
    end
end

for j = 1:length(BIons)
    matchIon = find(spectrum(:,1) > BIons(j)-obj.settings.MS2Tolerance & ...
        spectrum(:,1) < BIons(j)+obj.settings.MS2Tolerance);
    if ~isempty(matchIon)
       if numel(matchIon) > 1
          diff = abs(BIons(j)-spectrum(matchIon,1));
          min_diff = find(diff==min(diff));
          intensities = [intensities;spectrum(matchIon(min_diff),2)];  
       else
          intensities = [intensities;spectrum(matchIon,2)];
       end
    end
end

try
    intensities = sort(intensities,'descend');
    topN = sum(intensities(1:obj.settings.topNFragments,1));
catch
    try
        topN = sum(intensities(1:length(intensities)));
    catch
        topN = NaN;
    end
end

end

