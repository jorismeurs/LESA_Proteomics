function obj = plotMS2Data(obj,scanData)

if obj.settings.precursorRemoval == true
    charge = str2double(obj.output.peptideCharge(1));
    mz = calculateMZ(obj.output.peptideSequence,charge);
    maxDev = ppmDeviation(mz,obj.settings.MS1Tolerance);
    try
        idx = find(scanData(:,1) > mz-maxDev & scanData(:,1) < mz+maxDev);
        scanData(idx,:) = [];
    catch
        disp('Precursor not found');
    end
end

stem(scanData(:,1),scanData(:,2),'Color',repmat(0.5,1,3),'Marker','none');
xlim([min(scanData(:,1))-10 max(scanData(:,1))+10]);
ylim([0 max(scanData(:,2))+0.1*max(scanData(:,2))]); 

set(gcf,'Color','white');
set(gca,'FontName','Calibri','FontSize',14);
xlabel('\it{m/z}','interpreter','tex');
ylabel('Intensity');
title([obj.output.peptideSequence ' ' obj.output.peptideCharge])

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
              'Color',YColor,'Marker','none','LineWidth',obj.settings.width);
          text(scanData(matchIon(min_diff),1),scanData(matchIon(min_diff),2)+offset,...
              sprintf('y%d',j),'FontName','Calibri','FontSize',obj.settings.fontSize,...
              'HorizontalAlignment','Center',...
              'VerticalAlignment','Bottom',...
              'FontWeight','bold',...
              'FontAngle','italic');
       else
          offset = 20; 
          stem(scanData(matchIon,1),scanData(matchIon,2),...
              'Color',YColor,'Marker','none');
          text(scanData(matchIon,1),scanData(matchIon,2)+offset,...
              sprintf('y%d',j),'FontName','Calibri','FontSize',obj.settings.fontSize,...
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
              'Color',BColor,'Marker','none','LineWidth',obj.settings.width);
          text(scanData(matchIon(min_diff),1),scanData(matchIon(min_diff),2)+offset,...
              sprintf('b%d',j),'FontName','Calibri','FontSize',obj.settings.fontSize,...
              'HorizontalAlignment','Center',...
              'VerticalAlignment','Bottom',...
              'FontWeight','bold',...
              'FontAngle','italic');
       else
           offset = 0.1*scanData(matchIon,2);
          stem(scanData(matchIon,1),scanData(matchIon,2),...
              'Color',BColor,'Marker','none','LineWidth',obj.settings.width);
          text(scanData(matchIon,1),scanData(matchIon,2)+offset,...
              sprintf('b%d',j),'FontName','Calibri','FontSize',obj.settings.fontSize,...
              'HorizontalAlignment','Center',...
              'VerticalAlignment','Bottom',...
              'FontWeight','bold',...
              'FontAngle','italic');
       end
    end
end
hold off

if obj.settings.neutralLoss == true
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
                      'Color',YColor,'Marker','none','LineWidth',obj.settings.width);
                  text(scanData(matchIon(min_diff),1),scanData(matchIon(min_diff),2)+offset,...
                      sprintf('y%d-H_2O',j),'FontName','Calibri','FontSize',obj.settings.fontSize,...
                      'HorizontalAlignment','Center',...
                      'VerticalAlignment','Bottom',...
                      'FontWeight','bold',...
                      'FontAngle','italic');
               else
                  offset = 0.1*scanData(matchIon,2);
                  stem(scanData(matchIon,1),scanData(matchIon,2),...
                      'Color',YColor,'Marker','none','LineWidth',obj.settings.width);
                  text(scanData(matchIon,1),scanData(matchIon,2)+offset,...
                      sprintf('y%d-H_2O',j),'FontName','Calibri','FontSize',obj.settings.fontSize,...
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
                      'Color',YColor,'Marker','none','LineWidth',obj.settings.width);
                  text(scanData(matchIon(min_diff),1),scanData(matchIon(min_diff),2)+offset,...
                      sprintf('y%d-NH_3',j),'FontName','Calibri','FontSize',obj.settings.fontSize,...
                      'HorizontalAlignment','Center',...
                      'VerticalAlignment','Bottom',...
                      'FontWeight','bold',...
                      'FontAngle','italic');
               else
                  offset = 0.1*scanData(matchIon,2);
                  stem(scanData(matchIon,1),scanData(matchIon,2),...
                      'Color',YColor,'Marker','none','LineWidth',obj.settings.width);
                  text(scanData(matchIon,1),scanData(matchIon,2)+offset,...
                      sprintf('y%d-NH_3',j),'FontName','Calibri','FontSize',obj.settings.fontSize,...
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
                      'Color',YColor,'Marker','none','LineWidth',obj.settings.width);
                  text(scanData(matchIon(min_diff),1),scanData(matchIon(min_diff),2)+offset,...
                      sprintf('y%d-H_2O',j),'FontName','Calibri','FontSize',obj.settings.fontSize,...
                      'HorizontalAlignment','Center',...
                      'VerticalAlignment','Bottom',...
                      'FontWeight','bold',...
                      'FontAngle','italic');
               else
                  offset = 0.1*scanData(matchIon,2);
                  stem(scanData(matchIon,1),scanData(matchIon,2),...
                      'Color',YColor,'Marker','none','LineWidth',obj.settings.width);
                  text(scanData(matchIon,1),scanData(matchIon,2)+offset,...
                      sprintf('y%d-H_2O',j),'FontName','Calibri','FontSize',obj.settings.fontSize,...
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
                      'Color',YColor,'Marker','none','LineWidth',obj.settings.width);
                  text(scanData(matchIon(min_diff),1),scanData(matchIon(min_diff),2)+offset,...
                      sprintf('y%d-NH_3',j),'FontName','Calibri','FontSize',8,...
                      'HorizontalAlignment','Center',...
                      'VerticalAlignment','Bottom',...
                      'FontWeight','bold',...
                      'FontAngle','italic');
               else
                  offset = 0.1*scanData(matchIon,2);
                  stem(scanData(matchIon,1),scanData(matchIon,2),...
                      'Color',YColor,'Marker','none');
                  text(scanData(matchIon,1),scanData(matchIon,2)+offset,...
                      sprintf('y%d-NH_3',j),'FontName','Calibri','FontSize',8,...
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
end
set(gca,'box','off');
set(gcf,'Position',[50,50,1300,600]);

end

function mz = calculateMZ(sequence,charge)

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
for n = 1:numel(sequence)
   idx = find(ismember(aminoAcidMassList(:,1),sequence(n))==true);
   mz = [mz;aminoAcidMassList{idx,2}];
end
mz = sum(mz);

mz = mz+18.010565;

mz = (mz+(charge*1.007825)-(charge*0.00054858))/charge;


end