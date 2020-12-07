function plotMS1Data(obj)
clc
proteins = unique(obj.output.proteinList);
proteinSelection = questdlg('Single or all proteins?','Annotation',...
    'Single','All','Single');
switch proteinSelection
    case 'Single'
        for j = 1:length(proteins)
           fprintf('(%d) %s \n',j,proteins{j});
        end
        proteinIdx = input('Select protein: ');
        proteinLabel = input('Define label for protein: ','s');
        matchedRows = find(strcmp(proteins{proteinIdx},obj.output.proteinList));

        mzList = obj.output.mzList;
        tempList = [];
        for j = 1:length(mzList)
            tempMZ = str2num(cell2mat(obj.output.mzList(j)));
            tempList = [tempList;tempMZ];
        end
        mzList = tempList(matchedRows,1);
        charges = obj.output.chargeList(matchedRows,1);
        scanData = double(obj.output.scanData);
    case 'All'
        clc
        fprintf('PROVIDE LABELS FOR PLOT \n');
        for j = 1:length(proteins)
           proteinLabel{j} = input(sprintf('%s:  ',proteins{j}),'s'); 
        end
        scanData = double(obj.output.scanData);
    otherwise
        return
end

stem(scanData(:,1),scanData(:,2),'Color',repmat(0.5,1,3),'Marker','none');
if isempty(obj.settings.XLim)
    xlim([min(scanData(:,1))-10 max(scanData(:,1))+10]);
else
    xlim(obj.settings.XLim);
end
set(gcf,'Color','white');
set(gca,'FontName','Calibri','FontSize',14);
xlabel('\it{m/z}','interpreter','tex');
ylabel('Intensity');

precursorColor = [202/255,0/255,32/255];

switch proteinSelection
    case 'All'
        hold on
        for n = 1:length(proteins)
            for j = 1:length(obj.output.mzList)
                peakMatch = [];
                maxDev = ppmDeviation(str2num(obj.output.mzList{j}),obj.settings.MS1Tolerance);
                peakMatch = find(scanData(:,1) > str2num(obj.output.mzList{j})-maxDev & ...
                    scanData(:,1) < str2num(obj.output.mzList{j})+maxDev);
                if ~isempty(peakMatch)
                   if numel(peakMatch) > 1
                       diff = abs(scanData(peakMatch,1)-obj.output.mzList(j));
                       min_diff = find(diff==min(diff));
                       if scanData(peakMatch(min_diff),2)/max(scanData(:,2)) < 0.05
                           increaseVec = linspace(0.3,0.8,1000);
                           val = randperm(1000);
                           offset = increaseVec(val(1))*max(scanData(:,2));
                           stem(scanData(peakMatch(min_diff),1),scanData(peakMatch(min_diff),2),...
                               'Marker','none','Color',precursorColor,'LineWidth',obj.settings.width);
                           plot([scanData(peakMatch(min_diff),1) scanData(peakMatch(min_diff),1)],...
                               [scanData(peakMatch(min_diff),2) scanData(peakMatch(min_diff),2)+abs(offset)],...
                               '--','Color',repmat(0.1,1,3));
                           plot([scanData(peakMatch(min_diff),1) scanData(peakMatch(min_diff),1)+25],...
                               [scanData(peakMatch(min_diff),2)+abs(offset) scanData(peakMatch(min_diff),2)+abs(offset)],...
                               '--','Color',repmat(0.1,1,3));
                           text(scanData(peakMatch(min_diff),1)+26,scanData(peakMatch(min_diff),2)+abs(offset),...
                               proteinLabel{n},'FontSize',obj.settings.fontSize,'FontName','Calibri',...
                               'FontWeight','bold')
                       else
                           text(scanData(peakMatch(min_diff),1)+26,scanData(peakMatch(min_diff),2),...
                               proteinLabel{n},'FontSize',obj.settings.fontSize,'FontName','Calibri',...
                               'FontWeight','bold',...
                               'HorizontalAlignment','Center',...
                               'VerticalAlignment','Bottom')
                       end
                   else
                       stem(scanData(peakMatch,1),scanData(peakMatch,2),...
                           'Marker','none','Color',precursorColor,'LineWidth',obj.settings.width); 
                       if scanData(peakMatch,2)/max(scanData(:,2)) < 0.05
                           increaseVec = linspace(0.3,0.8,1000);
                           val = randperm(1000);
                           offset = increaseVec(val(1))*max(scanData(:,2));
                           plot([scanData(peakMatch,1) scanData(peakMatch,1)],...
                               [scanData(peakMatch,2) scanData(peakMatch,2)+abs(offset)],...
                               '--','Color',repmat(0.1,1,3));
                           plot([scanData(peakMatch,1) scanData(peakMatch,1)+25],...
                               [scanData(peakMatch,2)+abs(offset) scanData(peakMatch,2)+abs(offset)],...
                               '--','Color',repmat(0.1,1,3));
                           text(scanData(peakMatch,1)+26,scanData(peakMatch,2)+abs(offset),...
                               proteinLabel{n},'FontSize',obj.settings.fontSize,'FontName','Calibri',...
                               'FontWeight','bold')
                       else
                           text(scanData(peakMatch,1),scanData(peakMatch,2),...
                               proteinLabel{n},'FontSize',obj.settings.fontSize,'FontName','Calibri',...
                               'FontWeight','bold',...
                               'HorizontalAlignment','Center',...
                               'VerticalAlignment','Bottom');
                       end
                   end
                else
                   continue 
                end
            end
        end
    case 'Single'
        hold on
        obj.output.annotations = [];
        for j = 1:length(mzList)
            peakMatch = []; 
            maxDev = ppmDeviation(mzList(j),obj.settings.MS1Tolerance);
            peakMatch = find(scanData(:,1) > mzList(j)-maxDev & ...
                scanData(:,1) < mzList(j)+maxDev);
            if ~isempty(peakMatch)
               if numel(peakMatch) > 1
                   diff = abs(scanData(peakMatch,1)-mzList(j));
                   min_diff = find(diff==min(diff));
                   obj.output.annotations = [obj.output.annotations;scanData(peakMatch(min_diff),1),scanData(peakMatch(min_diff),2)];
                   fprintf('m/z %.4f (%.2e) \n',scanData(peakMatch(min_diff),1),scanData(peakMatch(min_diff),2));
                   if scanData(peakMatch(min_diff),2)/max(scanData(:,2)) < 0.05
                       increaseVec = linspace(0.3,0.8,1000);
                       val = randperm(1000);
                       offset = increaseVec(val(1))*max(scanData(:,2));
                       stem(scanData(peakMatch(min_diff),1),scanData(peakMatch(min_diff),2),...
                           'Marker','none','Color',precursorColor,'LineWidth',obj.settings.width);
                       plot([scanData(peakMatch(min_diff),1) scanData(peakMatch(min_diff),1)],...
                           [scanData(peakMatch(min_diff),2) scanData(peakMatch(min_diff),2)+abs(offset)],...
                           '--','Color',repmat(0.1,1,3));
                       plot([scanData(peakMatch(min_diff),1) scanData(peakMatch(min_diff),1)+25],...
                           [scanData(peakMatch(min_diff),2)+abs(offset) scanData(peakMatch(min_diff),2)+abs(offset)],...
                           '--','Color',repmat(0.1,1,3));
                       text(scanData(peakMatch(min_diff),1)+26,scanData(peakMatch(min_diff),2)+abs(offset),...
                           [proteinLabel ' (' charges{j} ')'],'FontSize',obj.settings.fontSize,'FontName','Calibri',...
                           'FontWeight','bold')
                   else
                       text(scanData(peakMatch(min_diff),1)+26,scanData(peakMatch(min_diff),2),...
                           [proteinLabel ' (' charges{j} ')'],'FontSize',obj.settings.fontSize,'FontName','Calibri',...
                           'FontWeight','bold',...
                           'HorizontalAlignment','Center',...
                           'VerticalAlignment','Bottom')
                   end
               else
                   obj.output.annotations = [obj.output.annotations;scanData(peakMatch,1),scanData(peakMatch,2)];
                   fprintf('m/z %.4f (%.2e) \n',scanData(peakMatch,1),scanData(peakMatch,2));
                   stem(scanData(peakMatch,1),scanData(peakMatch,2),...
                       'Marker','none','Color',precursorColor,'LineWidth',obj.settings.width); 
                   if scanData(peakMatch,2)/max(scanData(:,2)) < 0.05
                       increaseVec = linspace(0.3,0.8,1000);
                       val = randperm(1000);
                       offset = increaseVec(val(1))*max(scanData(:,2));
                       plot([scanData(peakMatch,1) scanData(peakMatch,1)],...
                           [scanData(peakMatch,2) scanData(peakMatch,2)+abs(offset)],...
                           '--','Color',repmat(0.1,1,3));
                       plot([scanData(peakMatch,1) scanData(peakMatch,1)+25],...
                           [scanData(peakMatch,2)+abs(offset) scanData(peakMatch,2)+abs(offset)],...
                           '--','Color',repmat(0.1,1,3));
                       text(scanData(peakMatch,1)+26,scanData(peakMatch,2)+abs(offset),...
                           [proteinLabel ' (' charges{j} ')'],'FontSize',obj.settings.fontSize,'FontName','Calibri',...
                           'FontWeight','bold')
                   else
                       text(scanData(peakMatch,1),scanData(peakMatch,2),...
                           [proteinLabel ' (' charges{j} ')'],'FontSize',obj.settings.fontSize,'FontName','Calibri',...
                           'FontWeight','bold',...
                           'HorizontalAlignment','Center',...
                           'VerticalAlignment','Bottom');
                   end
               end
            else
               continue 
            end
        end
        hold off
end
set(gca,'box','off');
set(gcf,'Position',[50,50,1300,600]);

end