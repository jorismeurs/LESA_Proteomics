function [yseries,bseries] = fragmentSequence(sequence)
% Generating theoretical m/z values for input peptide sequence
% Assumed charge state: 1+
% Included ions: b, y

% y-series
yseries = struct('y',[],'mz',[],'waterLoss',[],'ammoniaLoss',[]);
reversedSequence = fliplr(sequence);
for j = 1:numel(reversedSequence)
    tempFragment = reversedSequence(1:j);
    mz = calculateMZ(tempFragment,'y');
    yseries(j).y = sprintf('y%d',j);
    yseries(j).mz = mz;
end

% y-series + neutral losses
for j = 1:numel(reversedSequence)
    tempFragment = reversedSequence(1:j);
    
    % Calculate potential water loss (https://www.ru.nl/publish/pages/580560/denovopeptidesequencingtutorial.pdf)
    if contains(tempFragment,'E') || contains(tempFragment,'S') || contains(tempFragment,'T')
       mz = calculateMZ(tempFragment,'y');
       mz = mz-18.010565;
       yseries(j).waterLoss = mz; 
    else
       yseries(j).waterLoss = NaN; 
    end
    
     % Calculate potential ammonia loss (https://www.ru.nl/publish/pages/580560/denovopeptidesequencingtutorial.pdf)
    if contains(tempFragment,'K') || contains(tempFragment,'R') || contains(tempFragment,'Q') || contains(tempFragment,'N')
        mz = calculateMZ(tempFragment,'y');
        mz = mz-17.026548;
        yseries(j).ammoniaLoss = mz;
    else
        yseries(j).ammoniaLoss = NaN;
    end
end

% b-series
bseries = struct('b',[],'mz',[]);
for j = 1:numel(sequence)
    tempFragment = sequence(1:j);
    mz = calculateMZ(tempFragment,'b');
    bseries(j).b = sprintf('b%d',j);
    bseries(j).mz = mz;
end

% b-series + neutral losses
for j = 1:numel(sequence)
    tempFragment = sequence(1:j);
    
    % Calculate potential water loss (https://www.ru.nl/publish/pages/580560/denovopeptidesequencingtutorial.pdf)
    if contains(tempFragment,'E') || contains(tempFragment,'S') || contains(tempFragment,'T')
       mz = calculateMZ(tempFragment,'b');
       mz = mz-18.010565;
       bseries(j).waterLoss = mz; 
    else
       bseries(j).waterLoss = NaN; 
    end
    
     % Calculate potential ammonia loss (https://www.ru.nl/publish/pages/580560/denovopeptidesequencingtutorial.pdf)
    if contains(tempFragment,'K') || contains(tempFragment,'R') || contains(tempFragment,'Q') || contains(tempFragment,'N')
        mz = calculateMZ(tempFragment,'b');
        mz = mz-17.026548;
        bseries(j).ammoniaLoss = mz;
    else
        bseries(j).ammoniaLoss = NaN;
    end
end

end

function mz = calculateMZ(aminoAcids,ionseries)

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
for n = 1:numel(aminoAcids)
   idx = find(ismember(aminoAcidMassList(:,1),aminoAcids(n))==true);
   mz = [mz;aminoAcidMassList{idx,2}];
end
mz = sum(mz);

if isequal(ionseries,'y')
    mz = mz+18.010565+1.007825;
end
if isequal(ionseries,'b')
    mz = mz+1.007825;
end

end