function obj = plotLibraryData(obj,refSpec,sampleSpec,sequence,charge,corr)

sampleSpec(:,2) = (sampleSpec(:,2)./max(sampleSpec(:,2))).*100;
refSpec(:,2) = (refSpec(:,2)./max(refSpec(:,2))).*100;

if obj.settings.precursorRemoval == true
   precursorIon = getMZ(sequence,charge);
   tolerance = obj.settings.MS1Tolerance;
   maxDev = ppmDeviation(precursorIon,tolerance);
   precursorIndex = find(sampleSpec(:,1) > precursorIon-maxDev & ...
       sampleSpec(:,1) < precursorIon+maxDev);
   if ~isempty(precursorIndex)
      sampleSpec(precursorIndex,:) = []; 
   end
end

stem(sampleSpec(:,1),sampleSpec(:,2),'Color',repmat(0.7,1,3),'Marker','none');
xlim([min(sampleSpec(:,1))-10 max(sampleSpec(:,1))+10]);
ylim([-100 100]);
set(gcf,'Color','white');
set(gca,'FontName','Calibri','FontSize',14);
xlabel('\it{m/z}','interpreter','tex');
ylabel('Relative intensity (%)');
title([sprintf('%s %s ',sequence,charge) '(cos \theta: ' sprintf('%.4f',corr) ')'])

YIons = [obj.output.yIons.mz];
YColor = [202/255,0/255,32/255];
hold on
stem(refSpec(:,1),-refSpec(:,2),'Color',repmat(0.7,1,3),'Marker','none')
for j = 1:length(YIons)
    matchIon = find(sampleSpec(:,1) > YIons(j)-obj.settings.MS2Tolerance & ...
        sampleSpec(:,1) < YIons(j)+obj.settings.MS2Tolerance);
    if ~isempty(matchIon)
       if numel(matchIon) > 1
          diff = abs(YIons(j)-sampleSpec(matchIon,1));
          min_diff = find(diff==min(diff));
          offset = 5;
          stem(sampleSpec(matchIon(min_diff),1),sampleSpec(matchIon(min_diff),2),...
              'Color',YColor,'Marker','none');
          text(sampleSpec(matchIon(min_diff),1),sampleSpec(matchIon(min_diff),2)+offset,...
              sprintf('y%d',j),'FontName','Calibri','FontSize',8,...
              'HorizontalAlignment','Center',...
              'VerticalAlignment','Bottom',...
              'FontWeight','bold',...
              'FontAngle','italic');
       else
          offset = 5; 
          stem(sampleSpec(matchIon,1),sampleSpec(matchIon,2),...
              'Color',YColor,'Marker','none');
          text(sampleSpec(matchIon,1),sampleSpec(matchIon,2)+offset,...
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
    matchIon = find(sampleSpec(:,1) > BIons(j)-obj.settings.MS2Tolerance & ...
        sampleSpec(:,1) < BIons(j)+obj.settings.MS2Tolerance);
    if ~isempty(matchIon)
       if numel(matchIon) > 1
          diff = abs(BIons(j)-sampleSpec(matchIon,1));
          min_diff = find(diff==min(diff));
          offset = 5;
          stem(sampleSpec(matchIon(min_diff),1),sampleSpec(matchIon(min_diff),2),...
              'Color',YColor,'Marker','none');
          text(sampleSpec(matchIon(min_diff),1),sampleSpec(matchIon(min_diff),2)+offset,...
              sprintf('b%d',j),'FontName','Calibri','FontSize',8,...
              'HorizontalAlignment','Center',...
              'VerticalAlignment','Bottom',...
              'FontWeight','bold',...
              'FontAngle','italic');
       else
          offset = 5;
          stem(sampleSpec(matchIon,1),sampleSpec(matchIon,2),...
              'Color',YColor,'Marker','none');
          text(sampleSpec(matchIon,1),sampleSpec(matchIon,2)+offset,...
              sprintf('b%d',j),'FontName','Calibri','FontSize',8,...
              'HorizontalAlignment','Center',...
              'VerticalAlignment','Bottom',...
              'FontWeight','bold',...
              'FontAngle','italic');
          %obj.output.fragmentIonIntensity = [obj.output.fragmentIonIntensity;sampleSpec(matchIon,2)];
       end
    end
end
hold off

%set(gcf,'Position',[50,50,1300,600]);
set(gcf,'Position',[200,80,600,600]);
ytick = get(gca,'YTickLabel');
newLabels = [];
for j = 1:length(ytick)
    newLabels = [newLabels;abs(str2num(ytick{j}))];
end
newLabels = num2str(newLabels);
set(gca,'YTickLabel',newLabels);
end

