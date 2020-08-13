function R = cosineCorrelation(reference,sample)

% As described by Tabb, D. L., MacCoss, M. J., Wu, C. C., Anderson, S. D., & Yates, J. R. (2003). Similarity among Tandem Mass Spectra from Proteomic Experiments:  Detection, Significance, and Utility. Analytical Chemistry, 75(10), 2470–2477. doi:10.1021/ac026424o 

R = (sum(reference.*sample))/sqrt(sumsqr(reference)*sumsqr(sample));

end

