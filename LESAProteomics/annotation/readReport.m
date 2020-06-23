function reportData = readReport(obj)

fileID = fopen([obj.folder.identification '\A__B__1_Extended_PSM_Report.txt'],'r');
fileData = textscan(fileID,repmat('%s',1,25),'delimiter','\t','EmptyValue',Inf);
fclose(fileID);
fclose('all');

reportData = [];
for j = 1:size(fileData,2)
    reportData = [reportData,cellstr(fileData{j})];
end

end