dataDirs = {...
    '/home/jack/PhaseFeedforward/PhShftNetAnalyser_20150818/PFFPhaseShifterMeasurements_20150818_Shft1',...
    '/home/jack/PhaseFeedforward/PhShftNetAnalyser_20150818/PFFPhaseShifterMeasurements_20150818_Shft2',...
    '/home/jack/PhaseFeedforward/PhShftNetAnalyser_20150818/PFFPhaseShifterMeasurements_20150818_Shft3'...
};

dataSetNames = {'Shifter 1', 'Shifter 2', 'Shifter 3'};
dataSetColours = [0 0 1; 1 0 0; 0 0.5 0];
freqToUse = 3e9;

savePlots = 1;
saveDir = '/home/jack/PhaseFeedforward/PhShftNetAnalyser_20150818/Plots';

%% load data
addpath('../');
addpath('../../');

nDatsSets = length(dataDirs);

shifterSetting = cell(1,nDatsSets);
S11_Mag = cell(1,nDatsSets); % port 1 reflected power (dB)
S12_Mag = cell(1,nDatsSets); % port 1->2 insertion losses (dB)
S21_Mag = cell(1,nDatsSets); % port 2->1 insertion losses (dB)
S12_Phase = cell(1,nDatsSets); % port 1->2 phase (3GHz degrees)
S22_Mag = cell(1,nDatsSets); % port 2 reflected power (dB)
S21_Phase = cell(1,nDatsSets); % port 2->1 phase (3GHz degrees)

for i=1:nDatsSets
    
    fileList = dir([dataDirs{1} '/*.csv']);
    fileList = {fileList.name};
    
    fileValues = str2double(strrep(fileList,'.csv','')); % Assume file names of format value.csv. This strips .csv from the file names and converts the remaining string to a double.
    [shifterSetting{i},sortIndices] = sort(fileValues); % sort shifter settings in numerical order
    fileList = fileList(sortIndices); % sort file names in numerical order (using results of sort above)
    nFiles = length(fileList);
    
    S11_Mag{i} = NaN(1,nFiles); % port 1 reflected power (dB)
    S12_Mag{i} = NaN(1,nFiles); % port 1->2 insertion losses (dB)
    S21_Mag{i} = NaN(1,nFiles); % port 2->1 insertion losses (dB)
    S12_Phase{i} = NaN(1,nFiles); % port 1->2 phase (3GHz degrees)
    S22_Mag{i} = NaN(1,nFiles); % port 2 reflected power (dB)
    S21_Phase{i} = NaN(1,nFiles); % port 2->1 phase (3GHz degrees)
    
    for f=1:nFiles
        tmpData = importdata([dataDirs{i} '/' fileList{f}]);
        rowToExtract = find(tmpData.data(:,1)==freqToUse); % finds the row in first column of data (frquency values) that matches freqToUse
        
        % order of columns in files is assumed to be the same
        S11_Mag{i}(f) = tmpData.data(rowToExtract,2);
        S12_Mag{i}(f) = tmpData.data(rowToExtract,3);
        S21_Mag{i}(f) = tmpData.data(rowToExtract,4);
        S12_Phase{i}(f) = tmpData.data(rowToExtract,5);
        S22_Mag{i}(f) = tmpData.data(rowToExtract,6);
        S21_Phase{i}(f) = tmpData.data(rowToExtract,7);

    end
    
end

%% remove discontinuities in phase (i.e. let it go outside +/-180 degrees)
for i=1:nDatsSets
    phaseSteps = find(diff(S12_Phase{i})>0);        
    if (length(phaseSteps)>1)
        S12_Phase{i}(1:(phaseSteps(1))) = S12_Phase{i}(1:(phaseSteps(1)))+ 360;
        S12_Phase{i}((phaseSteps(2)+1):end) = S12_Phase{i}((phaseSteps(2)+1):end)- 360;
    else
        S12_Phase{i}((phaseSteps(1)+1):end) = S12_Phase{i}((phaseSteps(1)+1):end)- 360;   
    end
    
    phaseSteps = find(diff(S21_Phase{i})>0);
    if (length(phaseSteps)>1)
        S21_Phase{i}(1:(phaseSteps(1))) = S21_Phase{i}(1:(phaseSteps(1)))+ 360;
        S21_Phase{i}((phaseSteps(2)+1):end) = S21_Phase{i}((phaseSteps(2)+1):end)- 360;
    else
       S21_Phase{i}((phaseSteps(1)+1):end) = S21_Phase{i}((phaseSteps(1)+1):end)- 360; 
    end
end

%% fit phase vs. shifter setting

fitPhaseS21VsShifter = NaN(nDatsSets,2);
fitPhaseS12VsShifter = NaN(nDatsSets,2);
diffFitPhaseS21VsShifter = cell(1,nDatsSets);
diffFitPhaseS12VsShifter = cell(1,nDatsSets);

for i=1:nDatsSets
    fitPhaseS21VsShifter(i,:) = polyfit(shifterSetting{i}, S21_Phase{i},1);
    fitPhaseS12VsShifter(i,:) = polyfit(shifterSetting{i}, S12_Phase{i},1);
    
    fitY = fitPhaseS21VsShifter(i,1).*shifterSetting{i} + fitPhaseS21VsShifter(i,2);
    diffFitPhaseS21VsShifter{i} = S21_Phase{i}-fitY;
    
    fitY = fitPhaseS12VsShifter(i,1).*shifterSetting{i} + fitPhaseS12VsShifter(i,2);
    diffFitPhaseS12VsShifter{i} = S12_Phase{i}-fitY;

end

%% plots
format_plots;
close all;

figure;
for i=1:nDatsSets
    plot(shifterSetting{i},S11_Mag{i},'LineWidth',lineWidthBig);
    hold all;
end
title('S11 Amplitude')
xlabel('Shifter Setting [10x4GHz degrees]');
ylabel('Reflected Power [dB]')
legend(dataSetNames);
xlim([0 9000])
ylim([-60 -10]);
format_plots;
if savePlots; savePlot(saveDir,'S11_Mag'); end;

figure;
for i=1:nDatsSets
    plot(shifterSetting{i},S12_Mag{i},'LineWidth',lineWidthBig);
    hold all;
end
title('S12 Amplitude')
xlabel('Shifter Setting [10x4GHz degrees]');
ylabel('Insertion Losses [dB]')
legend(dataSetNames);
xlim([0 9000])
ylim([-1.5 -0.7])
format_plots;
if savePlots; savePlot(saveDir,'S12_Mag'); end;

figure;
for i=1:nDatsSets
    plot(shifterSetting{i},S21_Mag{i},'LineWidth',lineWidthBig);
    hold all;
end
title('S21 Amplitude')
xlabel('Shifter Setting [10x4GHz degrees]');
ylabel('Insertion Losses [dB]')
legend(dataSetNames);
xlim([0 9000])
ylim([-1.5 -0.7])
format_plots;
if savePlots; savePlot(saveDir,'S21_Mag'); end;

figure;
for i=1:nDatsSets
    plot(shifterSetting{i},S12_Phase{i},'LineWidth',lineWidthBig);
    hold all;
end
title('S12 Phase')
xlabel('Shifter Setting [10x4GHz degrees]');
ylabel('Phase [3GHz Degrees]')
legend(dataSetNames);
xlim([0 9000])
%ylim([-180 180]);
format_plots;
if savePlots; savePlot(saveDir,'S12_Phase'); end;

figure;
for i=1:nDatsSets
    plot(shifterSetting{i},diffFitPhaseS12VsShifter{i},'LineWidth',lineWidthBig);
    hold all;
end
title('S12 Phase: Difference to Linear Fit')
xlabel('Shifter Setting [10x4GHz degrees]');
ylabel('Phase Difference [3GHz Degrees]')
legend(dataSetNames);
xlim([0 9000])
format_plots;
if savePlots; savePlot(saveDir,'S12_DiffPhase'); end;

figure;
for i=1:nDatsSets
    plot(shifterSetting{i},S22_Mag{i},'LineWidth',lineWidthBig);
    hold all;
end
title('S22 Amplitude')
xlabel('Shifter Setting [10x4GHz degrees]');
ylabel('Reflected Power [dB]')
legend(dataSetNames);
xlim([0 9000])
ylim([-60 -10]);
format_plots;
if savePlots; savePlot(saveDir,'S22_Mag'); end;

figure;
for i=1:nDatsSets
    plot(shifterSetting{i},S21_Phase{i},'LineWidth',lineWidthBig);
    hold all;
end
title('S21 Phase')
xlabel('Shifter Setting [10x4GHz degrees]');
ylabel('Phase [3GHz Degrees]')
legend(dataSetNames);
xlim([0 9000])
%ylim([-180 180]);
format_plots;
if savePlots; savePlot(saveDir,'S21_Phase'); end;

figure;
for i=1:nDatsSets
    plot(shifterSetting{i},diffFitPhaseS21VsShifter{i},'LineWidth',lineWidthBig);
    hold all;
end
title('S21 Phase: Difference to Linear Fit')
xlabel('Shifter Setting [10x4GHz degrees]');
ylabel('Phase Difference [3GHz Degrees]')
legend(dataSetNames);
xlim([0 9000])
format_plots;
if savePlots; savePlot(saveDir,'S21_DiffPhase'); end;

