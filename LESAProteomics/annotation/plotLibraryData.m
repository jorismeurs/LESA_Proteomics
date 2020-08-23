function obj = plotLibraryData(obj,refSpec,sampleSpec,originalSpec,fragments,sequence,charge,corr)

if obj.settings.precursorRemoval == true
   precursorIon = getMZ(sequence,charge);
   tolerance = obj.settings.MS1Tolerance;
   maxDev = ppmDeviation(precursorIon,tolerance);
   precursorIndex = find(originalSpec(:,1) > precursorIon-maxDev & ...
       originalSpec(:,1) < precursorIon+maxDev);
   if ~isempty(precursorIndex)
      originalSpec(precursorIndex,:) = []; 
   end
end

originalSpec(:,2) = (originalSpec(:,2)./max(originalSpec(:,2))).*100;
refSpec(:,2) = (refSpec(:,2)./max(refSpec(:,2))).*100;
sampleSpec(:,2) = (sampleSpec(:,2)./max(sampleSpec(:,2))).*100;

stem(originalSpec(:,1),originalSpec(:,2),'Color',repmat(0.7,1,3),'Marker','none');
xlim([min(originalSpec(:,1))-10 max(originalSpec(:,1))+10]);
ylim([-100 100]);
set(gcf,'Color','white');
set(gca,'FontName','Calibri','FontSize',14);
xlabel('\it{m/z}','interpreter','tex');
ylabel('Relative intensity (%)');
title([sprintf('%s %s ',sequence,charge) '(cos \theta: ' sprintf('%.4f',corr) ')'])

YColor = [202/255,0/255,32/255];
BColor = [5/255,113/255,176/255];
hold on
stem(refSpec(:,1),-refSpec(:,2),'Color',repmat(0.7,1,3),'Marker','none')
for j = 1:length(sampleSpec)
   if sampleSpec(j,2) > 0 
       tempFragment = fragments(j,1);
       if contains(tempFragment,'b')
           stem(sampleSpec(j,1),sampleSpec(j,2),'Marker','none','Color',BColor);
           text(sampleSpec(j,1),sampleSpec(j,2),tempFragment,...
               'HorizontalAlignment','Center','VerticalAlignment','Bottom',...
               'FontAngle','italic',...
               'FontWeight','bold',...
               'FontName','Calibri');
       elseif contains(tempFragment,'y')
           stem(sampleSpec(j,1),sampleSpec(j,2),'Marker','none','Color',YColor);
           text(sampleSpec(j,1),sampleSpec(j,2),tempFragment,...
               'HorizontalAlignment','Center','VerticalAlignment','Bottom',...
               'FontAngle','italic',...
               'FontWeight','bold',...
               'FontName','Calibri');
       end
   else
       continue
   end
end
hold off


set(gcf,'Position',[50,50,1300,600]);
%set(gcf,'Position',[200,80,600,600]);
ytick = get(gca,'YTickLabel');
newLabels = [];
for j = 1:length(ytick)
    newLabels = [newLabels;abs(str2num(ytick{j}))];
end
newLabels = num2str(newLabels);
set(gca,'YTickLabel',newLabels);
end

