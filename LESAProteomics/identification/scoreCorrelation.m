function ID = scoreCorrelation(MGFStruct,libraryData)

ID = struct('protein',[],...
    'sequence',[],...
    'sample',[],...
    'reference',[],...
    'R',[]);

count = 0;
tolerance = 20; % ppm
for j = 1:length(MGFStruct.scan)
    scanData = MGFStruct.scan(j).scanData;
    for n = 1:length(libraryData)
       librarySpectrum = libraryData(n).spectrum;
       libraryMZ = librarySpectrum(:,1);
       intensityMatrix = zeros(length(libraryMZ),2);
       sampleMZ = [];
       for k = 1:length(libraryMZ)           
           maxDev = ppmDeviation(libraryMZ(k),tolerance); 
           idx = find(scanData(:,1) > libraryMZ(k)-maxDev & ...
               scanData(:,1) < libraryMZ(k)+maxDev);
           if ~isempty(idx)
               if numel(idx) > 1
                   diff = libraryMZ(k)-scanData(idx,1);
                   min_diff = find(diff==min(diff));
                   sampleMZ = [sampleMZ;scanData(idx(min_diff),1)];
                   intensityMatrix(k,1) = scanData(idx(min_diff),2);
               else
                   sampleMZ = [sampleMZ;scanData(idx,1)];
                   intensityMatrix(k,1) = scanData(idx,2);
               end
           else
               sampleMZ = [sampleMZ;libraryMZ(k)];
           end
       end
       intensityMatrix(:,1) = (intensityMatrix(:,1).*100)./max(intensityMatrix(:,1));
       intensityMatrix(:,2) = (librarySpectrum(:,2).*100)./max(librarySpectrum(:,2));
       referenceSpectrum = [libraryMZ,intensityMatrix(:,2)];
       sampleSpectrum = [sampleMZ,intensityMatrix(:,1)];
       if ~isempty(sampleSpectrum) && ~isempty(referenceSpectrum)
           R = corr(referenceSpectrum,sampleSpectrum);
           if R >= 0.95
                count = count+1;
                ID(count).protein = libraryData(n).protein;
                ID(count).sequence = libraryData(n).sequence;
                ID(count).sample = sampleSpectrum;
                ID(count).reference = referenceSpectrum;
                ID(count).R = R;
           else
               continue
           end
       end
    end
end

end
