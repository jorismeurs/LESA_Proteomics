function jitterplot(data,y,varargin)
%
% INTRODUCTION
% Function to create jitter plots showing the mean, mean - 1 standard
% deviation, mean + 1 standard deviation and the distribution of the
% individual data points
%
% INPUT
% jitterplot(data,y)
% data:     N x 1 vector containing the data points to be plotted
% y:        N x 1 vector representing the corresponding groups of the data
%           points
% 
% jitterplot(...,'meanWidth',meanWidthValue)
% Value for the width of the line to indicate the mean in the plot. Default
% is 0.2
%
% jitterplot(...,'stdWidth',stdWithValue)
% Value for the width of the lines to indicate the standard deviation in
% the plot. Default is 0.3
%
% jitterplot(...,'scatterSize',scatterSizeValue)
% Value for the marker size to show the individual data points. Default
% value is 20
%
% Created:      8/7/2020
% Version:      1.1
% Last update:  9/7/2020
% Developer:    Joris Meurs

if nargin < 2
   error('Not enough input arguments'); 
end

meanWidth = [];
stdWidth = [];
scatterSize = [];

for j = 1:length(varargin)
   if isequal(varargin{j},'meanWidth')
      meanWidth = varargin{j+1}; 
   end
   if isequal(varargin{j},'stdWidth')
      stdWidth = varargin{j+1}; 
   end
   if isequal(varargin{j},'scatterSize')
      scatterSize = varargin{j+1}; 
   end
end

if isempty(meanWidth)
   meanWidth = 0.2; 
end

if isempty(stdWidth)
   stdWidth = 0.3; 
end

if isempty(scatterSize)
   scatterSize = 20; 
end

if length(data) ~= length(y)
   error('Data and group vector should be the same length'); 
end

uniqueGroups = unique(y,'stable');

% Supported colors
% (https://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12)
c = [
    228,26,28
    55,126,184
    77,175,74
    152,78,163
    ];
c = c./255;

hold on
for j = 1:length(uniqueGroups)
    % Calculate mean and standard deviation for group
    mu = nanmean(data(uniqueGroups(j)==y));
    sd = nanstd(data(uniqueGroups(j)==y));

    
    % Plot mean line
    plot([j-meanWidth/2 j+meanWidth/2],[mu mu],'-','Color',c(j,:),'LineWidth',2);
    
    % Plot mean-standard deviation
    plot([j j],[mu-sd mu],'-','Color',c(j,:),'LineWidth',2);
    
    % Plot mean+standard deviation
    plot([j j],[mu mu+sd],'-','Color',c(j,:),'LineWidth',2);
    
    % Plot mean-standard deviation end
    plot([j-stdWidth/2 j+stdWidth/2],[mu-sd mu-sd],'-','Color',c(j,:),'LineWidth',2);
    
    % Plot mean+standard deviation end
    plot([j-stdWidth/2 j+stdWidth/2],[mu+sd mu+sd],'-','Color',c(j,:),'LineWidth',2);
    
    % Scatter plot for data points
    x = 0.05.*rand(length(find(uniqueGroups(j)==y)),1)+j;
    scatter(x,data(uniqueGroups(j)==y),scatterSize,'filled',...
        'MarkerFaceColor','k');
end
hold off

end