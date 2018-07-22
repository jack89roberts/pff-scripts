bpmPrefixes = cell(1,nBPMs);
for i=1:nBPMs
    bpmPrefixes{i} = bpmNames{i}(1:2);
end
firstCL = find(ismember(bpmPrefixes,'CL')); 
try
    firstCL = firstCL(1);
catch
    firstCL = -1;
end
firstCT = find(ismember(bpmPrefixes,'CT')); 
try
    firstCT = firstCT(1);
catch
    firstCT = -1;
end
firstCR = find(ismember(bpmPrefixes,'CR')); 
try
    firstCR = firstCR(1);
catch
    firstCR = -1;
end
firstCC = find(ismember(bpmPrefixes,'CC')); 
try
    firstCC = firstCC(1);
catch
    firstCC = -1;
end
firstCB = find(ismember(bpmPrefixes,'CB')); 
try
    firstCB = firstCB(1);
catch
    firstCB = -1;
end

xMax = ceil(nBPMs/10)*10;

set(gcf,'PaperPositionMode','auto'); % forces figures to be saved with same aspect ratio
set(gcf,'units','normalized','position',[0 0 0.8 0.5])
set(gca,'XTick',0:5:xMax)
set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)
hold all;

if (firstCT ~= -1)
    plot([firstCT firstCT], [yMin yMax],'k','LineWidth',2)
end
if (firstCR ~= -1)
    plot([firstCR firstCR], [yMin yMax],'k','LineWidth',2)
end
if (firstCC ~=-1)
    plot([firstCC firstCC], [yMin yMax],'k','LineWidth',2)
end
if (firstCB ~=-1)
    plot([firstCB firstCB], [yMin yMax],'k','LineWidth',2)
end
plot([0 xMax], [corrMon2Mon3 corrMon2Mon3],'k','LineWidth',2);
plot([0 xMax], [-corrMon2Mon3 -corrMon2Mon3],'k','LineWidth',2);
grid on;
xlabel('BPM Index')
ylabel('Correlation')
ylim([yMin yMax])

% Section text box labels
% annotation(gcf,'textbox',...
%     [0.271633237822348 0.163434525791914 0.0422129468188442 0.0903273789699944],...
%     'String',{'CT'},...
%     'FontSize',18,...
%     'FontName','Helvetica',...
%     'EdgeColor',[1 1 1]);
% 
% annotation(gcf,'textbox',...
%     [0.468767908309452 0.163434525791915 0.0435260575288962 0.0903273789699944],...
%     'String',{'CR'},...
%     'FontSize',18,...
%     'FontName','Helvetica',...
%     'EdgeColor',[1 1 1]);
% 
% annotation(gcf,'textbox',...
%     [0.610888252148995 0.163434525791914 0.0435794980810494 0.0903273789699944],...
%     'String',{'CC'},...
%     'LineStyle','none',...
%     'FontSize',18,...
%     'FontName','Helvetica',...
%     'EdgeColor',[1 1 1]);
% 
% annotation(gcf,'textbox',...
%     [0.791977077363893 0.163434525791914 0.0433886389662163 0.0903273789699944],...
%     'String',{'CB'},...
%     'FontSize',18,...
%     'FontName','Helvetica',...
%     'EdgeColor',[1 1 1]);
% 
% annotation(gcf,'textbox',...
%     [0.135243553008595 0.163434525791912 0.0413731667135785 0.0903273789699944],...
%     'String',{'CL'},...
%     'LineStyle','none',...
%     'FontSize',18,...
%     'FontName','Helvetica',...
%     'EdgeColor',[1 1 1]);

% Create textbox
strSections = 'Sections Included: ';
if (firstCL ~= -1)
    strSections = [strSections 'CL '];
end
if (firstCT ~= -1)
    strSections = [strSections 'CT '];
end
if (firstCR ~= -1)
    strSections = [strSections 'CR '];
end
if (firstCC ~=-1)
    strSections = [strSections 'CC '];
end
if (firstCB ~=-1)
    strSections = [strSections 'CB '];
end

annotation(gcf,'textbox',...
    [0.589843750000002 0.165805557405783 0.304425039988382 0.084305553705328],...
    'String', strSections,...
    'FontSize',18,...
    'FontName','Helvetica',...
    'FitBoxToText','off',...
    'BackgroundColor',[1 1 1]);
