function obj = plotMS2Data(obj,scanData)

stem(scanData(:,1),scanData(:,2),'Color',repmat(0.7,1,3),'Marker','none');
xlim([100 1500]);
set(gcf,'Color','white');
set(gca,'FontName','Calibri');
xlabel('\it{m/z}','interpreter','tex');
ylabel('Intensity');

YIons = [obj.output.yIons.mz];
YColor = [202/255,0/255,32/255];
hold on
for j = 1:length(YIons)
    matchIon = find(scanData(:,1) > YIons(j)-obj.settings.MS2Tolerance & ...
        scanData(:,1) < YIons(j)+obj.settings.MS2Tolerance);
    if ~isempty(matchIon)
       if numel(matchIon) > 1
          diff = YIons(j)-scanData(matchIon,1);
          min_diff = find(diff==min(diff));
          offset = 0.1*scanData(matchIon(min_diff),2);
          stem(scanData(matchIon(min_diff),1),scanData(matchIon(min_diff),2),...
              'Color',YColor,'Marker','none');
          text(scanData(matchIon(min_diff),1),scanData(matchIon(min_diff),2)+offset,...
              sprintf('y%d',j),'FontName','Calibri','FontSize',8,...
              'HorizontalAlignment','Center',...
              'VerticalAlignment','Bottom',...
              'FontWeight','bold',...
              'FontAngle','italic');
       else
           offset = 0.1*scanData(matchIon,2);
          stem(scanData(matchIon,1),scanData(matchIon,2),...
              'Color',YColor,'Marker','none');
          text(scanData(matchIon,1),scanData(matchIon,2)+offset,...
              sprintf('y%d',j),'FontName','Calibri','FontSize',8,...
              'HorizontalAlignment','Center',...
              'VerticalAlignment','Bottom',...
              'FontWeight','bold',...
              'FontAngle','italic');
       end
    end
end
hold off

BIons = [obj.output.bIons.mz];
YColor = [5/255,113/255,176/255];
hold on
for j = 1:length(BIons)
    matchIon = find(scanData(:,1) > BIons(j)-obj.settings.MS2Tolerance & ...
        scanData(:,1) < BIons(j)+obj.settings.MS2Tolerance);
    if ~isempty(matchIon)
       if numel(matchIon) > 1
          diff = BIons(j)-scanData(matchIon,1);
          min_diff = find(diff==min(diff));
          offset = 0.1*scanData(matchIon(min_diff),2);
          stem(scanData(matchIon(min_diff),1),scanData(matchIon(min_diff),2),...
              'Color',YColor,'Marker','none');
          text(scanData(matchIon(min_diff),1),scanData(matchIon(min_diff),2)+offset,...
              sprintf('b%d',j),'FontName','Calibri','FontSize',8,...
              'HorizontalAlignment','Center',...
              'VerticalAlignment','Bottom',...
              'FontWeight','bold',...
              'FontAngle','italic');
       else
           offset = 0.1*scanData(matchIon,2);
          stem(scanData(matchIon,1),scanData(matchIon,2),...
              'Color',YColor,'Marker','none');
          text(scanData(matchIon,1),scanData(matchIon,2)+offset,...
              sprintf('b%d',j),'FontName','Calibri','FontSize',8,...
              'HorizontalAlignment','Center',...
              'VerticalAlignment','Bottom',...
              'FontWeight','bold',...
              'FontAngle','italic');
       end
    end
end
hold off
set(gcf,'Position',[100,100,1000,600]);