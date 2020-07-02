function obj = plotMS2Data(obj,scanData)

stem(scanData(:,1),scanData(:,2),'Color',repmat(0.7,1,3),'Marker','none');
xlim([min(scanData(:,1))-10 max(scanData(:,1))+10]);
%ylim([-100 100]);
set(gcf,'Color','white');
set(gca,'FontName','Calibri','FontSize',14);
xlabel('\it{m/z}','interpreter','tex');
ylabel('Intensity');
title([obj.output.peptideSequence ' ' obj.output.precursorCharge])

YIons = [obj.output.yIons.mz];
YColor = [202/255,0/255,32/255];
hold on
for j = 1:length(YIons)
    matchIon = find(scanData(:,1) > YIons(j)-obj.settings.MS2Tolerance & ...
        scanData(:,1) < YIons(j)+obj.settings.MS2Tolerance);
    if ~isempty(matchIon)
       if numel(matchIon) > 1
          diff = abs(YIons(j)-scanData(matchIon,1));
          min_diff = find(diff==min(diff));
          offset = 20;
          stem(scanData(matchIon(min_diff),1),scanData(matchIon(min_diff),2),...
              'Color',YColor,'Marker','none');
          text(scanData(matchIon(min_diff),1),scanData(matchIon(min_diff),2)+offset,...
              sprintf('y%d',j),'FontName','Calibri','FontSize',8,...
              'HorizontalAlignment','Center',...
              'VerticalAlignment','Bottom',...
              'FontWeight','bold',...
              'FontAngle','italic');
       else
          offset = 20; 
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
BColor = [5/255,113/255,176/255];
hold on
for j = 1:length(BIons)
    matchIon = find(scanData(:,1) > BIons(j)-obj.settings.MS2Tolerance & ...
        scanData(:,1) < BIons(j)+obj.settings.MS2Tolerance);
    if ~isempty(matchIon)
       if numel(matchIon) > 1
          diff = abs(BIons(j)-scanData(matchIon,1));
          min_diff = find(diff==min(diff));
          offset = 0.1*scanData(matchIon(min_diff),2);
          stem(scanData(matchIon(min_diff),1),scanData(matchIon(min_diff),2),...
              'Color',BColor,'Marker','none');
          text(scanData(matchIon(min_diff),1),scanData(matchIon(min_diff),2)+offset,...
              sprintf('b%d',j),'FontName','Calibri','FontSize',8,...
              'HorizontalAlignment','Center',...
              'VerticalAlignment','Bottom',...
              'FontWeight','bold',...
              'FontAngle','italic');
       else
           offset = 0.1*scanData(matchIon,2);
          stem(scanData(matchIon,1),scanData(matchIon,2),...
              'Color',BColor,'Marker','none');
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

% Neutral losses y-ions
YLoss = [obj.output.yIons.waterLoss];
hold on
for j = 1:length(YLoss)
   if isnan(YLoss(j))
       continue
   else
       matchIon = find(scanData(:,1) > YLoss(j)-obj.settings.MS2Tolerance & ...
           scanData(:,1) < YLoss(j)+obj.settings.MS2Tolerance);
       if ~isempty(matchIon)
           if numel(matchIon) > 1
              diff = abs(YLoss(j)-scanData(matchIon,1));
              min_diff = find(diff==min(diff));
              offset = 0.1*scanData(matchIon(min_diff),2);
              stem(scanData(matchIon(min_diff),1),scanData(matchIon(min_diff),2),...
                  'Color',YColor,'Marker','none');
              text(scanData(matchIon(min_diff),1),scanData(matchIon(min_diff),2)+offset,...
                  sprintf('y%d-H_2O',j),'FontName','Calibri','FontSize',8,...
                  'HorizontalAlignment','Center',...
                  'VerticalAlignment','Bottom',...
                  'FontWeight','bold',...
                  'FontAngle','italic');
           else
              offset = 0.1*scanData(matchIon,2);
              stem(scanData(matchIon,1),scanData(matchIon,2),...
                  'Color',YColor,'Marker','none');
              text(scanData(matchIon,1),scanData(matchIon,2)+offset,...
                  sprintf('y%d-H_2O',j),'FontName','Calibri','FontSize',8,...
                  'HorizontalAlignment','Center',...
                  'VerticalAlignment','Bottom',...
                  'FontWeight','bold',...
                  'FontAngle','italic');
           end
       else
           continue
       end
   end
end
hold off

YLoss = [obj.output.yIons.ammoniaLoss];
hold on
for j = 1:length(YLoss)
   if isnan(YLoss(j))
       continue
   else
       matchIon = find(scanData(:,1) > YLoss(j)-obj.settings.MS2Tolerance & ...
           scanData(:,1) < YLoss(j)+obj.settings.MS2Tolerance);
       if ~isempty(matchIon)
           if numel(matchIon) > 1
              diff = abs(YLoss(j)-scanData(matchIon,1));
              min_diff = find(diff==min(diff));
              offset = 0.1*scanData(matchIon(min_diff),2);
              stem(scanData(matchIon(min_diff),1),scanData(matchIon(min_diff),2),...
                  'Color',YColor,'Marker','none');
              text(scanData(matchIon(min_diff),1),scanData(matchIon(min_diff),2)+offset,...
                  sprintf('y%d-H_2O',j),'FontName','Calibri','FontSize',8,...
                  'HorizontalAlignment','Center',...
                  'VerticalAlignment','Bottom',...
                  'FontWeight','bold',...
                  'FontAngle','italic');
           else
              offset = 0.1*scanData(matchIon,2);
              stem(scanData(matchIon,1),scanData(matchIon,2),...
                  'Color',YColor,'Marker','none');
              text(scanData(matchIon,1),scanData(matchIon,2)+offset,...
                  sprintf('y%d-H_2O',j),'FontName','Calibri','FontSize',8,...
                  'HorizontalAlignment','Center',...
                  'VerticalAlignment','Bottom',...
                  'FontWeight','bold',...
                  'FontAngle','italic');
           end
       else
           continue
       end
   end
end
hold off

% Neutral losses b-ions
BLoss = [obj.output.bIons.waterLoss];
hold on
for j = 1:length(BLoss)
   if isnan(BLoss(j))
       continue
   else
       matchIon = find(scanData(:,1) > BLoss(j)-obj.settings.MS2Tolerance & ...
           scanData(:,1) < BLoss(j)+obj.settings.MS2Tolerance);
       if ~isempty(matchIon)
           if numel(matchIon) > 1
              diff = abs(BLoss(j)-scanData(matchIon,1));
              min_diff = find(diff==min(diff));
              offset = 0.1*scanData(matchIon(min_diff),2);
              stem(scanData(matchIon(min_diff),1),scanData(matchIon(min_diff),2),...
                  'Color',YColor,'Marker','none');
              text(scanData(matchIon(min_diff),1),scanData(matchIon(min_diff),2)+offset,...
                  sprintf('y%d-H_2O',j),'FontName','Calibri','FontSize',8,...
                  'HorizontalAlignment','Center',...
                  'VerticalAlignment','Bottom',...
                  'FontWeight','bold',...
                  'FontAngle','italic');
           else
              offset = 0.1*scanData(matchIon,2);
              stem(scanData(matchIon,1),scanData(matchIon,2),...
                  'Color',YColor,'Marker','none');
              text(scanData(matchIon,1),scanData(matchIon,2)+offset,...
                  sprintf('y%d-H_2O',j),'FontName','Calibri','FontSize',8,...
                  'HorizontalAlignment','Center',...
                  'VerticalAlignment','Bottom',...
                  'FontWeight','bold',...
                  'FontAngle','italic');
           end
       else
           continue
       end
   end
end
hold off

BLoss = [obj.output.bIons.ammoniaLoss];
hold on
for j = 1:length(BLoss)
   if isnan(BLoss(j))
       continue
   else
       matchIon = find(scanData(:,1) > BLoss(j)-obj.settings.MS2Tolerance & ...
           scanData(:,1) < BLoss(j)+obj.settings.MS2Tolerance);
       if ~isempty(matchIon)
           if numel(matchIon) > 1
              diff = abs(BLoss(j)-scanData(matchIon,1));
              min_diff = find(diff==min(diff));
              offset = 0.1*scanData(matchIon(min_diff),2);
              stem(scanData(matchIon(min_diff),1),scanData(matchIon(min_diff),2),...
                  'Color',YColor,'Marker','none');
              text(scanData(matchIon(min_diff),1),scanData(matchIon(min_diff),2)+offset,...
                  sprintf('y%d-H_2O',j),'FontName','Calibri','FontSize',8,...
                  'HorizontalAlignment','Center',...
                  'VerticalAlignment','Bottom',...
                  'FontWeight','bold',...
                  'FontAngle','italic');
           else
              offset = 0.1*scanData(matchIon,2);
              stem(scanData(matchIon,1),scanData(matchIon,2),...
                  'Color',YColor,'Marker','none');
              text(scanData(matchIon,1),scanData(matchIon,2)+offset,...
                  sprintf('y%d-H_2O',j),'FontName','Calibri','FontSize',8,...
                  'HorizontalAlignment','Center',...
                  'VerticalAlignment','Bottom',...
                  'FontWeight','bold',...
                  'FontAngle','italic');
           end
       else
           continue
       end
   end
end
hold off

set(gcf,'Position',[50,50,1300,600]);