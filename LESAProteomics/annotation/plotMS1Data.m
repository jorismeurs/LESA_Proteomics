function plotMS1Data(obj)

scanData = double(obj.output.scanData);

stem(scanData(:,1),scanData(:,2),'Color',repmat(0.7,1,3),'Marker','none');
xlim([min(scanData(:,1))-10 max(scanData(:,1))+10]);
set(gcf,'Color','white');
set(gca,'FontName','Calibri','FontSize',14);
xlabel('\it{m/z}','interpreter','tex');
ylabel('Intensity');

precursorColor = [202/255,0/255,32/255];

hold on
for j = 1:length(obj.output.mzList)
    maxDev = ppmDeviation(str2num(obj.output.mzList{j}),obj.settings.MS1Tolerance);
    peakMatch = find(scanData(:,1) > str2num(obj.output.mzList{j})-maxDev & ...
        scanData(:,1) < str2num(obj.output.mzList{j})+maxDev);
    if ~isempty(peakMatch)
       if numel(peakMatch) > 1
           diff = scanData(peakMatch,1)-obj.output.mzList(j);
           min_diff = find(diff==min(diff));
           if scanData(peakMatch(min_diff),2)/max(scanData(:,2)) < 0.05
               increaseVec = linspace(0.3,0.8,1000);
               val = randperm(1000);
               offset = increaseVec(val(1))*max(scanData(:,2));
               stem(scanData(peakMatch(min_diff),1),scanData(peakMatch(min_diff),2),...
                   'Marker','none','Color',precursorColor);
               plot([scanData(peakMatch(min_diff),1) scanData(peakMatch(min_diff),1)],...
                   [scanData(peakMatch(min_diff),2) scanData(peakMatch(min_diff),2)+abs(offset)],...
                   'k:');
               plot([scanData(peakMatch(min_diff),1) scanData(peakMatch(min_diff),1)+25],...
                   [scanData(peakMatch(min_diff),2)+abs(offset) scanData(peakMatch(min_diff),2)+abs(offset)],...
                   'k:');
               text(scanData(peakMatch(min_diff),1)+26,scanData(peakMatch(min_diff),2)+abs(offset),...
                   obj.output.proteinList{j},'FontSize',9,'FontName','Calibri',...
                   'FontWeight','bold')
           else
               text(scanData(peakMatch(min_diff),1)+26,scanData(peakMatch(min_diff),2),...
                   obj.output.proteinList{j},'FontSize',9,'FontName','Calibri',...
                   'FontWeight','bold')
           end
       else
           stem(scanData(peakMatch,1),scanData(peakMatch,2),...
               'Marker','none','Color',precursorColor); 
           if scanData(peakMatch,2)/max(scanData(:,2)) < 0.05
               increaseVec = linspace(0.3,0.8,1000);
               val = randperm(1000);
               offset = increaseVec(val(1))*max(scanData(:,2));
               plot([scanData(peakMatch,1) scanData(peakMatch,1)],...
                   [scanData(peakMatch,2) scanData(peakMatch,2)+abs(offset)],...
                   'k:');
               plot([scanData(peakMatch,1) scanData(peakMatch,1)+25],...
                   [scanData(peakMatch,2)+abs(offset) scanData(peakMatch,2)+abs(offset)],...
                   'k:');
               text(scanData(peakMatch,1)+26,scanData(peakMatch,2)+abs(offset),...
                   obj.output.proteinList{j},'FontSize',9,'FontName','Calibri',...
                   'FontWeight','bold')
           else
               text(scanData(peakMatch,1)+26,scanData(peakMatch,2),...
                   obj.output.proteinList{j},'FontSize',9,'Calibri',...
                   'FontWeight','bold')
           end
       end
    end
end
set(gcf,'Position',[50,50,1300,600]);

end