function R = cosineCorrelation(reference,sample)

% Original developer: R.G. Bettinardi
% Ruggero G. Bettinardi (2020). getCosineSimilarity(x,y) (https://www.mathworks.com/matlabcentral/fileexchange/62978-getcosinesimilarity-x-y), MATLAB Central File Exchange. Retrieved July 1, 2020.

xy   = dot(reference,sample);
nx   = norm(reference);
ny   = norm(sample);
nxny = nx*ny;
R   = xy/nxny;

end

