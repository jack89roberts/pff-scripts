%% plots

%%%%%%%%%% Scatter plots for max correlation (H,S,V) acros BPMs %%%%%%%%%%
% Mon3 - BPMH
if (displayAll)
    figure;
else
    hold off;
end
[~,i] = max(abs(corrMon3_BPMH));
scatter(meanPulsePhase(3,:),meanBPMH{i});
title(sprintf('Mon3-%sH: %.2f',bpmNames{i},corrMon3_BPMH(i)) ) 
xlabel('Mon3 Phase [degrees]');
ylabel('BPM Position [mm]');
if (savePlots)
    saveStr = sprintf('%s/corr_Mon3_BPMH_%s',saveDir,bpmNames{i});
    print([saveStr '.png'],'-dpng');
    savefig([saveStr '.fig']);
end

% Mon3 - BPMV
if (useVertical)
    [~,i] = max(abs(corrMon3_BPMV));
    isNaNBPM = isnan(meanBPMV(i,:));
    isGoodPulse = ~(isNaNBPM | isNaN3');
    goodBPM = meanBPMV(i,isGoodPulse);
    goodMon3 = meanMon3Phase(isGoodPulse);
    tmpCorr = corrcoef(goodMon3,goodBPM);
    corrMon3_BPMV(i) = tmpCorr(1,2);
    if (displayAll)
        figure;
    else
        hold off;
    end
    scatter(goodMon3,goodBPM);
    title(sprintf('Mon3-%sV: %.2f',bpmNames{i},corrMon3_BPMV(i)) )
    xlabel('Mon3 Phase [degrees]');
    ylabel('BPM Position [mm]');
    if (savePlots)
        saveStr = sprintf('%s/corr_Mon3_BPMV_%s',saveDir,bpmNames{i});
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end
end

% Mon3 - BPMS
[~,i] = max(abs(corrMon3_BPMS));
isNaNBPM = isnan(meanBPMS(i,:));
isGoodPulse = ~(isNaNBPM | isNaN3');
goodBPM = meanBPMS(i,isGoodPulse);
goodMon3 = meanMon3Phase(isGoodPulse);
tmpCorr = corrcoef(goodMon3,goodBPM);
corrMon3_BPMS(i) = tmpCorr(1,2);
if (displayAll)
    figure;
else
    hold off;
end
scatter(goodMon3,goodBPM);
title(sprintf('Mon3-%sS: %.2f',bpmNames{i},corrMon3_BPMS(i)) )
xlabel('Mon3 Phase [degrees]');
ylabel('BPM Transmission [V]');
if (savePlots)
    saveStr = sprintf('%s/corr_Mon3_BPMS_%s',saveDir,bpmNames{i});
    print([saveStr '.png'],'-dpng');
    savefig([saveStr '.fig']);
end

% Mon2 - BPMH
[~,i] = max(abs(corrMon2_BPMH));
isNaNBPM = isnan(meanBPMH(i,:));
isGoodPulse = ~(isNaNBPM | isNaN2');
goodBPM = meanBPMH(i,isGoodPulse);
goodMon2 = meanMon2Phase(isGoodPulse);
tmpCorr = corrcoef(goodMon2,goodBPM);
corrMon2_BPMH(i) = tmpCorr(1,2);
if (displayAll)
    figure;
else
    hold off;
end
scatter(goodMon2,goodBPM);
title(sprintf('Mon2-%sH: %.2f',bpmNames{i},corrMon2_BPMH(i)) )
xlabel('Mon2 Phase [degrees]');
ylabel('BPM Position [mm]');
if (savePlots)
    saveStr = sprintf('%s/corr_Mon2_BPMH_%s',saveDir,bpmNames{i});
    print([saveStr '.png'],'-dpng');
    savefig([saveStr '.fig']);
end

% Mon2 - BPMV
if (useVertical)
    [~,i] = max(abs(corrMon2_BPMV));
    isNaNBPM = isnan(meanBPMV(i,:));
    isGoodPulse = ~(isNaNBPM | isNaN2');
    goodBPM = meanBPMV(i,isGoodPulse);
    goodMon2 = meanMon2Phase(isGoodPulse);
    tmpCorr = corrcoef(goodMon2,goodBPM);
    corrMon2_BPMV(i) = tmpCorr(1,2);
    if (displayAll)
        figure;
    else
        hold off;
    end
    scatter(goodMon2,goodBPM);
    title(sprintf('Mon2-%sV: %.2f',bpmNames{i},corrMon2_BPMV(i)) )
    xlabel('Mon2 Phase [degrees]');
    ylabel('BPM Position [mm]');
    if (savePlots)
        saveStr = sprintf('%s/corr_Mon2_BPMV_%s',saveDir,bpmNames{i});
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end
end

% Mon2 - BPMS
[~,i] = max(abs(corrMon2_BPMS));
isNaNBPM = isnan(meanBPMS(i,:));
isGoodPulse = ~(isNaNBPM | isNaN2');
goodBPM = meanBPMS(i,isGoodPulse);
goodMon2 = meanMon2Phase(isGoodPulse);
tmpCorr = corrcoef(goodMon2,goodBPM);
corrMon2_BPMS(i) = tmpCorr(1,2);
if (displayAll)
    figure;
else
    hold off;
end
scatter(goodMon2,goodBPM);
title(sprintf('Mon2-%sS: %.2f',bpmNames{i},corrMon2_BPMS(i)) )
xlabel('Mon2 Phase [degrees]');
ylabel('BPM Transmission [V]');
if (savePlots)
    saveStr = sprintf('%s/corr_Mon2_BPMS_%s',saveDir,bpmNames{i});
    print([saveStr '.png'],'-dpng');
    savefig([saveStr '.fig']);
end

%%%%%%%%%% Correlation vs. BPM %%%%%%%%%%

yMin = -1;
yMax = 1;

if (useVertical)    
    if (displayAll)
        figure;
    else
        hold off;
    end
    plot(corrMon1_BPMV,'LineWidth',2);    
    hold all;
    plot(corrMon2_BPMV,'LineWidth',2);
    plot(corrMon3_BPMV,'LineWidth',2);
    title('VERTICAL POSITION: Correlation with Phase')
    legend('Mon2 (CT)','Mon3 (CB)');
    format_plots_corrVsBPM;
    if (savePlots)
        saveStr = sprintf('%s/corrBPMV_Phase',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end
end

if (displayAll)
    figure;
else
    hold off;
end
plot(corrMon2_BPMH,'b','LineWidth',2);
hold all;
plot(corrMon3_BPMH,'r','LineWidth',2);
title('HORIZONTAL POSITION: Correlation with Phase')
legend('Mon2 (CT)','Mon3 (CB)');xlabel('BPM Index')
format_plots_corrVsBPM;
if (savePlots)
    saveStr = sprintf('%s//corrBPMH_Phase',saveDir);
    print([saveStr '.png'],'-dpng');
    savefig([saveStr '.fig']);
end

if (displayAll)
    figure;
else
    hold off;
end
plot(corrMon2_BPMS,'b','LineWidth',2);
hold all;
plot(corrMon3_BPMS,'r','LineWidth',2);
title('BPM TRANSMISSION: Correlation with Phase')
legend('Mon2 (CT)','Mon3 (CB)');
format_plots_corrVsBPM;
if (savePlots)
    saveStr = sprintf('%s//corrBPMS_Phase',saveDir);
    print([saveStr '.png'],'-dpng');
    savefig([saveStr '.fig']);
end

%%%%%% Correlations: BPM CL502/CT285 signals and VERTICAL POSITION %%%%%%
if (useVertical)
    green = [0, 0.8, 0];
    yMin = -1;
    yMax = 1;

    if (useCL502)    
        if (displayAll)
            figure;
        else
            hold off;
        end
        plot(corrMon2_BPMV,'b','LineWidth',2);
        hold all;
        plot(corrMon3_BPMV,'r','LineWidth',2);
        plot(corrBPMCL502H_BPMV, 'Color',green,'LineWidth',2)
        plot([bpmCL502Index bpmCL502Index],[yMin yMax],'--','Color',green,'LineWidth',2);
        title('Correlation with Vertical Position in all BPMs')
        legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CL502H');
        format_plots_corrVsBPM;
        if (savePlots)
            saveStr = sprintf('%s//corrBPMV_CL502H',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end
        
        if (displayAll)
            figure;
        else
            hold off;
        end
        plot(corrMon2_BPMV,'b','LineWidth',2);
        hold all;
        plot(corrMon3_BPMV,'r','LineWidth',2);
        plot(corrBPMCL502S_BPMV, 'Color',green,'LineWidth',2)
        plot([bpmCL502Index bpmCL502Index],[yMin yMax],'--','Color',green,'LineWidth',2);
        title('Correlation with Vertical Position in all BPMs')
        legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CL502S');
        format_plots_corrVsBPM;
        if (savePlots)
            saveStr = sprintf('%s//corrBPMV_CL502S',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        if (displayAll)
            figure;
        else
            hold off;
        end
        plot(corrMon2_BPMV,'b','LineWidth',2);
        hold all;
        plot(corrMon3_BPMV,'r','LineWidth',2);
        plot(corrBPMCL502V_BPMV, 'Color',green,'LineWidth',2)
        plot([bpmCL502Index bpmCL502Index],[yMin yMax],'--','Color',green,'LineWidth',2);
        title('Correlation with Vertical Position in all BPMs')
        legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CL502V');
        format_plots_corrVsBPM;
        if (savePlots)
            saveStr = sprintf('%s//corrBPMV_CL502V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end
    end
    
    if (useCT285)
        if (displayAll)
            figure;
        else
            hold off;
        end
        plot(corrMon2_BPMV,'b','LineWidth',2);
        hold all;
        plot(corrMon3_BPMV,'r','LineWidth',2);
        plot(corrBPMCT285H_BPMV, 'Color',green,'LineWidth',2)
        plot([bpmCT285Index bpmCT285Index],[yMin yMax],'--','Color',green,'LineWidth',2);
        title('Correlation with Vertical Position in all BPMs')
        grid on;
        legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CT285H');
        format_plots_corrVsBPM;
        if (savePlots)
            saveStr = sprintf('%s//corrBPMV_CT285H',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        if (displayAll)
            figure;
        else
            hold off;
        end
        plot(corrMon2_BPMV,'b','LineWidth',2);
        hold all;
        plot(corrMon3_BPMV,'r','LineWidth',2);
        plot(corrBPMCT285S_BPMV, 'Color',green,'LineWidth',2)
        plot([bpmCT285Index bpmCT285Index],[yMin yMax],'--','Color',green,'LineWidth',2);
        title('Correlation with Vertical Position in all BPMs')
        legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CT285S');
        format_plots_corrVsBPM;
        if (savePlots)
            saveStr = sprintf('%s//corrBPMV_CT285S',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        if (displayAll)
            figure;
        else
            hold off;
        end
        plot(corrMon2_BPMV,'b','LineWidth',2);
        hold all;
        plot(corrMon3_BPMV,'r','LineWidth',2);
        plot(corrBPMCT285V_BPMV, 'Color',green,'LineWidth',2)
        plot([bpmCT285Index bpmCT285Index],[yMin yMax],'--','Color',green,'LineWidth',2);
        title('Correlation with Vertical Position in all BPMs')
        legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CT285V');
        format_plots_corrVsBPM;
        if (savePlots)
            saveStr = sprintf('%s//corrBPMV_CT285V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end 
    end
end

%%%%%% Correlations: BPM CL502/CT285 signals and HORIZONTAL POSITION %%%%%%
green = [0, 0.8, 0];
yMin = -1;
yMax = 1;

if (useCL502)    
    if (displayAll)
        figure;
    else
        hold off;
    end
    plot(corrMon2_BPMH,'b','LineWidth',2);
    hold all;
    plot(corrMon3_BPMH,'r','LineWidth',2);
    plot(corrBPMCL502H_BPMH, 'Color',green,'LineWidth',2)
    plot([bpmCL502Index bpmCL502Index],[yMin yMax],'--','Color',green,'LineWidth',2);
    title('Correlation with Horizontal Position in all BPMs')
    legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CL502H');
    format_plots_corrVsBPM;
    if (savePlots)
        saveStr = sprintf('%s//corrBPMH_CL502H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    if (displayAll)
        figure;
    else
        hold off;
    end
    plot(corrMon2_BPMH,'b','LineWidth',2);
    hold all;
    plot(corrMon3_BPMH,'r','LineWidth',2);
    plot(corrBPMCL502S_BPMH, 'Color',green,'LineWidth',2)
    plot([bpmCL502Index bpmCL502Index],[yMin yMax],'--','Color',green,'LineWidth',2);
    title('Correlation with Horizontal Position in all BPMs')
    legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CL502S');
    format_plots_corrVsBPM;
    if (savePlots)
        saveStr = sprintf('%s//corrBPMH_CL502S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    if (useVertical)
        if (displayAll)
            figure;
        else
            hold off;
        end
        plot(corrMon2_BPMH,'b','LineWidth',2);
        hold all;
        plot(corrMon3_BPMH,'r','LineWidth',2);
        plot(corrBPMCL502V_BPMH, 'Color',green,'LineWidth',2)
        plot([bpmCL502Index bpmCL502Index],[yMin yMax],'--','Color',green,'LineWidth',2);
        title('Correlation with Horizontal Position in all BPMs')
        legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CL502V');
        format_plots_corrVsBPM;
        if (savePlots)
            saveStr = sprintf('%s//corrBPMH_CL502V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end
    end
end

if (useCT285)
    if (displayAll)
        figure;
    else
        hold off;
    end
    plot(corrMon2_BPMH,'b','LineWidth',2);
    hold all;
    plot(corrMon3_BPMH,'r','LineWidth',2);
    plot(corrBPMCT285H_BPMH, 'Color',green,'LineWidth',2)
    plot([bpmCT285Index bpmCT285Index],[yMin yMax],'--','Color',green,'LineWidth',2);
    title('Correlation with Horizontal Position in all BPMs')
    legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CT285H');
    format_plots_corrVsBPM;
    if (savePlots)
        saveStr = sprintf('%s//corrBPMH_CT285H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    if (displayAll)
        figure;
    else
        hold off;
    end
    plot(corrMon2_BPMH,'b','LineWidth',2);
    hold all;
    plot(corrMon3_BPMH,'r','LineWidth',2);
    plot(corrBPMCT285S_BPMH, 'Color',green,'LineWidth',2)
    plot([bpmCT285Index bpmCT285Index],[yMin yMax],'--','Color',green,'LineWidth',2);
    title('Correlation with Horizontal Position in all BPMs')
    legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CT285S');
    format_plots_corrVsBPM;
    if (savePlots)
        saveStr = sprintf('%s/corrBPMH_CT285S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    if (useVertical)
        if (displayAll)
            figure;
        else
            hold off;
        end
        plot(corrMon2_BPMH,'b','LineWidth',2);
        hold all;
        plot(corrMon3_BPMH,'r','LineWidth',2);
        plot(corrBPMCT285V_BPMH, 'Color',green,'LineWidth',2)
        plot([bpmCT285Index bpmCT285Index],[yMin yMax],'--','Color',green,'LineWidth',2);
        title('Correlation with Horizontal Position in all BPMs')
        legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CT285V');
        format_plots_corrVsBPM;
        if (savePlots)
            saveStr = sprintf('%s/corrBPMH_CT285V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end
    end
end

%%%%%% Correlations: BPM CL502/CT285 signals and TRANSMISSION %%%%%%
green = [0, 0.8, 0];
yMin = -1;
yMax = 1;

if (useCL502)
    if (displayAll)
        figure;
    else
        hold off;
    end
    plot(corrMon2_BPMS,'b','LineWidth',2);
    hold all;
    plot(corrMon3_BPMS,'r','LineWidth',2);
    plot(corrBPMCL502H_BPMS, 'Color',green,'LineWidth',2)
    plot([bpmCL502Index bpmCL502Index],[yMin yMax],'--','Color',green,'LineWidth',2);
    title('Correlation with Transmission in all BPMs')
    legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CL502H');
    format_plots_corrVsBPM;
    if (savePlots)
        saveStr = sprintf('%s/corrBPMS_CL502H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    if (displayAll)
        figure;
    else
        hold off;
    end
    plot(corrMon2_BPMS,'b','LineWidth',2);
    hold all;
    plot(corrMon3_BPMS,'r','LineWidth',2);
    plot(corrBPMCL502S_BPMS, 'Color',green,'LineWidth',2)
    plot([bpmCL502Index bpmCL502Index],[yMin yMax],'--','Color',green,'LineWidth',2);
    title('Correlation with Transmission in all BPMs')
    legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CL502S');
    format_plots_corrVsBPM;
    if (savePlots)
        saveStr = sprintf('%s/corrBPMS_CL502S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    if (useVertical)
        if (displayAll)
            figure;
        else
            hold off;
        end
        plot(corrMon2_BPMS,'b','LineWidth',2);
        hold all;
        plot(corrMon3_BPMS,'r','LineWidth',2);
        plot(corrBPMCL502V_BPMS, 'Color',green,'LineWidth',2)
        plot([bpmCL502Index bpmCL502Index],[yMin yMax],'--','Color',green,'LineWidth',2);
        title('Correlation with Transmission in all BPMs')
        legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CL502V');
        format_plots_corrVsBPM;
        if (savePlots)
            saveStr = sprintf('%s/corrBPMS_CL502V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end
    end
end

if (useCT285)
    if (displayAll)
        figure;
    else
        hold off;
    end
    plot(corrMon2_BPMS,'b','LineWidth',2);
    hold all;
    plot(corrMon3_BPMS,'r','LineWidth',2);
    plot(corrBPMCT285H_BPMS, 'Color',green,'LineWidth',2)
    plot([bpmCT285Index bpmCT285Index],[yMin yMax],'--','Color',green,'LineWidth',2);
    title('Correlation with Transmission in all BPMs')
    legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CT285H');
    format_plots_corrVsBPM;
    if (savePlots)
        saveStr = sprintf('%s/corrBPMS_CT285H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    if (displayAll)
        figure;
    else
        hold off;
    end
    plot(corrMon2_BPMS,'b','LineWidth',2);
    hold all;
    plot(corrMon3_BPMS,'r','LineWidth',2);
    plot(corrBPMCT285S_BPMS, 'Color',green,'LineWidth',2)
    plot([bpmCT285Index bpmCT285Index],[yMin yMax],'--','Color',green,'LineWidth',2);
    title('Correlation with Transmission in all BPMs')
    legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CT285S');
    format_plots_corrVsBPM;
    if (savePlots)
        saveStr = sprintf('%s/corrBPMS_CT285S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    if (useVertical)
        if (displayAll)
            figure;
        else
            hold off;
        end
        plot(corrMon2_BPMS,'b','LineWidth',2);
        hold all;
        plot(corrMon3_BPMS,'r','LineWidth',2);
        plot(corrBPMCT285V_BPMS, 'Color',green,'LineWidth',2)
        plot([bpmCT285Index bpmCT285Index],[yMin yMax],'--','Color',green,'LineWidth',2);
        title('Correlation with Transmission in all BPMs')
        legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CT285V');
        format_plots_corrVsBPM;
        if (savePlots)
            saveStr = sprintf('%s/corrBPMS_CT285V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end
    end
end

%% Correlation between phase correlation and bpm correlation

%% Mon3
if (displayAll)
    figure;
else
    hold off;
end
if (useRefBPM2)
    scatter(corrMon3_BPMH,corrRefBPM2H_BPMH);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon3 Correlation')
    ylabel('BPM CT 285H Correlation')
    title({'Correlation with BPMH: Mon3 vs. RefBPM2H', sprintf('(correl = %.2f)',bpmH_CORR_Mon3_RefBPM2H)});
    if (savePlots)
        saveStr = sprintf('%s/bpmH_CORR_Mon3_RefBPM2H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    scatter(corrMon3_BPMS,corrRefBPM2H_BPMS);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon3 Correlation')
    ylabel('BPM CT 285H Correlation')
    title({'Correlation with BPMS: Mon3 vs. RefBPM2H', sprintf('(correl = %.2f)',bpmS_CORR_Mon3_RefBPM2H)});
    if (savePlots)
        saveStr = sprintf('%s/bpmS_CORR_Mon3_RefBPM2H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    scatter(corrMon3_BPMH,corrRefBPM2S_BPMH);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon3 Correlation')
    ylabel('BPM CT 285S Correlation')
    title({'Correlation with BPMH: Mon3 vs. RefBPM2S', sprintf('(correl = %.2f)',bpmH_CORR_Mon3_RefBPM2S)});
    if (savePlots)
        saveStr = sprintf('%s/bpmH_CORR_Mon3_RefBPM2S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    scatter(corrMon3_BPMS,corrRefBPM2S_BPMS);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon3 Correlation')
    ylabel('BPM CT 285S Correlation')
    title({'Correlation with BPMS: Mon3 vs. RefBPM2S', sprintf('(correl = %.2f)',bpmS_CORR_Mon3_RefBPM2S)});
    if (savePlots)
        saveStr = sprintf('%s/bpmS_CORR_Mon3_RefBPM2S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end
end

if (useRefBPM1)
    scatter(corrMon3_BPMH,corrRefBPM1H_BPMH);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon3 Correlation')
    ylabel('BPM CL 502H Correlation')
    title({'Correlation with BPMH: Mon3 vs. RefBPM1H', sprintf('(correl = %.2f)',bpmH_CORR_Mon3_RefBPM1H)});
    if (savePlots)
        saveStr = sprintf('%s/bpmH_CORR_Mon3_RefBPM1H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    scatter(corrMon3_BPMS,corrRefBPM1H_BPMS);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon3 Correlation')
    ylabel('BPM CL 502H Correlation')
    title({'Correlation with BPMS: Mon3 vs. RefBPM1H', sprintf('(correl = %.2f)',bpmS_CORR_Mon3_RefBPM1H)});
    if (savePlots)
        saveStr = sprintf('%s/bpmS_CORR_Mon3_RefBPM1H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end


    scatter(corrMon3_BPMH,corrRefBPM1S_BPMH);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon3 Correlation')
    ylabel('BPM CL 502S Correlation')
    title({'Correlation with BPMH: Mon3 vs. RefBPM1S', sprintf('(correl = %.2f)',bpmH_CORR_Mon3_RefBPM1S)});
    if (savePlots)
        saveStr = sprintf('%s/bpmH_CORR_Mon3_RefBPM1S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    scatter(corrMon3_BPMS,corrRefBPM1S_BPMS);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon3 Correlation')
    ylabel('BPM CL 502S Correlation')
    title({'Correlation with BPMS: Mon3 vs. RefBPM1S', sprintf('(correl = %.2f)',bpmS_CORR_Mon3_RefBPM1S)});
    if (savePlots)
        saveStr = sprintf('%s/bpmS_CORR_Mon3_RefBPM1S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end
end

if (useVertical)
    if (useRefBPM2)
        scatter(corrMon3_BPMV,corrRefBPM2H_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon3 Correlation')
        ylabel('BPM CT 285H Correlation')
        title({'Correlation with BPMV: Mon3 vs. RefBPM2H', sprintf('(correl = %.2f)',bpmV_CORR_Mon3_RefBPM2H)});
        if (savePlots)
            saveStr = sprintf('%s/bpmV_CORR_Mon3_RefBPM2H',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon3_BPMH,corrRefBPM2V_BPMH);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon3 Correlation')
        ylabel('BPM CT 285V Correlation')
        title({'Correlation with BPMH: Mon3 vs. RefBPM2V', sprintf('(correl = %.2f)',bpmH_CORR_Mon3_RefBPM2V)});
        if (savePlots)
            saveStr = sprintf('%s/bpmH_CORR_Mon3_RefBPM2V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon3_BPMS,corrRefBPM2V_BPMS);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon3 Correlation')
        ylabel('BPM CT 285V Correlation')
        title({'Correlation with BPMS: Mon3 vs. RefBPM2V', sprintf('(correl = %.2f)',bpmS_CORR_Mon3_RefBPM2V)});
        if (savePlots)
            saveStr = sprintf('%s/bpmS_CORR_Mon3_RefBPM2V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon3_BPMV,corrRefBPM2V_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon3 Correlation')
        ylabel('BPM CT 285V Correlation')
        title({'Correlation with BPMV: Mon3 vs. RefBPM2V', sprintf('(correl = %.2f)',bpmV_CORR_Mon3_RefBPM2V)});
        if (savePlots)
            saveStr = sprintf('%s/bpmV_CORR_Mon3_RefBPM2V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon3_BPMV,corrRefBPM2S_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon3 Correlation')
        ylabel('BPM CT 285S Correlation')
        title({'Correlation with BPMV: Mon3 vs. RefBPM2S', sprintf('(correl = %.2f)',bpmV_CORR_Mon3_RefBPM2S)});
        if (savePlots)
            saveStr = sprintf('%s/bpmV_CORR_Mon3_RefBPM2S',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end
    end
    
    if (useRefBPM1)
        scatter(corrMon3_BPMV,corrRefBPM1H_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon3 Correlation')
        ylabel('BPM CL 502H Correlation')
        title({'Correlation with BPMV: Mon3 vs. RefBPM1H', sprintf('(correl = %.2f)',bpmV_CORR_Mon3_RefBPM1H)});
        if (savePlots)
            saveStr = sprintf('%s/bpmV_CORR_Mon3_RefBPM1H',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon3_BPMH,corrRefBPM1V_BPMH);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon3 Correlation')
        ylabel('BPM CL 502V Correlation')
        title({'Correlation with BPMH: Mon3 vs. RefBPM1V', sprintf('(correl = %.2f)',bpmH_CORR_Mon3_RefBPM1V)});
        if (savePlots)
            saveStr = sprintf('%s/bpmH_CORR_Mon3_RefBPM1V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon3_BPMS,corrRefBPM1V_BPMS);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon3 Correlation')
        ylabel('BPM CL 502V Correlation')
        title({'Correlation with BPMS: Mon3 vs. RefBPM1V', sprintf('(correl = %.2f)',bpmS_CORR_Mon3_RefBPM1V)});
        if (savePlots)
            saveStr = sprintf('%s/bpmS_CORR_Mon3_RefBPM1V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon3_BPMV,corrRefBPM1V_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon3 Correlation')
        ylabel('BPM CL 502V Correlation')
        title({'Correlation with BPMV: Mon3 vs. RefBPM1V', sprintf('(correl = %.2f)',bpmV_CORR_Mon3_RefBPM1V)});
        if (savePlots)
            saveStr = sprintf('%s/bpmV_CORR_Mon3_RefBPM1V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon3_BPMV,corrRefBPM1S_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon3 Correlation')
        ylabel('BPM CL 502S Correlation')
        title({'Correlation with BPMV: Mon3 vs. RefBPM1S', sprintf('(correl = %.2f)',bpmV_CORR_Mon3_RefBPM1S)});
        if (savePlots)
            saveStr = sprintf('%s/bpmV_CORR_Mon3_RefBPM1S',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end
    end

end

%% Mon2

if (useRefBPM2)
    scatter(corrMon2_BPMH,corrRefBPM2H_BPMH);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon2 Correlation')
    ylabel('BPM CT 285H Correlation')
    title({'Correlation with BPMH: Mon2 vs. RefBPM2H', sprintf('(correl = %.2f)',bpmH_CORR_Mon2_RefBPM2H)});
    if (savePlots)
        saveStr = sprintf('%s/bpmH_CORR_Mon2_RefBPM2H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    scatter(corrMon2_BPMS,corrRefBPM2H_BPMS);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon2 Correlation')
    ylabel('BPM CT 285H Correlation')
    title({'Correlation with BPMS: Mon2 vs. RefBPM2H', sprintf('(correl = %.2f)',bpmS_CORR_Mon2_RefBPM2H)});
    if (savePlots)
        saveStr = sprintf('%s/bpmS_CORR_Mon2_RefBPM2H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    scatter(corrMon2_BPMH,corrRefBPM2S_BPMH);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon2 Correlation')
    ylabel('BPM CT 285S Correlation')
    title({'Correlation with BPMH: Mon2 vs. RefBPM2S', sprintf('(correl = %.2f)',bpmH_CORR_Mon2_RefBPM2S)});
    if (savePlots)
        saveStr = sprintf('%s/bpmH_CORR_Mon2_RefBPM2S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    scatter(corrMon2_BPMS,corrRefBPM2S_BPMS);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon2 Correlation')
    ylabel('BPM CT 285S Correlation')
    title({'Correlation with BPMS: Mon2 vs. RefBPM2S', sprintf('(correl = %.2f)',bpmS_CORR_Mon2_RefBPM2S)});
    if (savePlots)
        saveStr = sprintf('%s/bpmS_CORR_Mon2_RefBPM2S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end
end

if (useRefBPM1)
    scatter(corrMon2_BPMH,corrRefBPM1H_BPMH);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon2 Correlation')
    ylabel('BPM CL 502H Correlation')
    title({'Correlation with BPMH: Mon2 vs. RefBPM1H', sprintf('(correl = %.2f)',bpmH_CORR_Mon2_RefBPM1H)});
    if (savePlots)
        saveStr = sprintf('%s/bpmH_CORR_Mon2_RefBPM1H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    scatter(corrMon2_BPMS,corrRefBPM1H_BPMS);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon2 Correlation')
    ylabel('BPM CL 502H Correlation')
    title({'Correlation with BPMS: Mon2 vs. RefBPM1H', sprintf('(correl = %.2f)',bpmS_CORR_Mon2_RefBPM1H)});
    if (savePlots)
        saveStr = sprintf('%s/bpmS_CORR_Mon2_RefBPM1H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end


    scatter(corrMon2_BPMH,corrRefBPM1S_BPMH);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon2 Correlation')
    ylabel('BPM CL 502S Correlation')
    title({'Correlation with BPMH: Mon2 vs. RefBPM1S', sprintf('(correl = %.2f)',bpmH_CORR_Mon2_RefBPM1S)});
    if (savePlots)
        saveStr = sprintf('%s/bpmH_CORR_Mon2_RefBPM1S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    scatter(corrMon2_BPMS,corrRefBPM1S_BPMS);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon2 Correlation')
    ylabel('BPM CL 502S Correlation')
    title({'Correlation with BPMS: Mon2 vs. RefBPM1S', sprintf('(correl = %.2f)',bpmS_CORR_Mon2_RefBPM1S)});
    if (savePlots)
        saveStr = sprintf('%s/bpmS_CORR_Mon2_RefBPM1S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end
end

if (useVertical)
    if (useRefBPM2)
        scatter(corrMon2_BPMV,corrRefBPM2H_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon2 Correlation')
        ylabel('BPM CT 285H Correlation')
        title({'Correlation with BPMV: Mon2 vs. RefBPM2H', sprintf('(correl = %.2f)',bpmV_CORR_Mon2_RefBPM2H)});
        if (savePlots)
            saveStr = sprintf('%s/bpmV_CORR_Mon2_RefBPM2H',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon2_BPMH,corrRefBPM2V_BPMH);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon2 Correlation')
        ylabel('BPM CT 285V Correlation')
        title({'Correlation with BPMH: Mon2 vs. RefBPM2V', sprintf('(correl = %.2f)',bpmH_CORR_Mon2_RefBPM2V)});
        if (savePlots)
            saveStr = sprintf('%s/bpmH_CORR_Mon2_RefBPM2V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon2_BPMS,corrRefBPM2V_BPMS);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon2 Correlation')
        ylabel('BPM CT 285V Correlation')
        title({'Correlation with BPMS: Mon2 vs. RefBPM2V', sprintf('(correl = %.2f)',bpmS_CORR_Mon2_RefBPM2V)});
        if (savePlots)
            saveStr = sprintf('%s/bpmS_CORR_Mon2_RefBPM2V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon2_BPMV,corrRefBPM2V_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon2 Correlation')
        ylabel('BPM CT 285V Correlation')
        title({'Correlation with BPMV: Mon2 vs. RefBPM2V', sprintf('(correl = %.2f)',bpmV_CORR_Mon2_RefBPM2V)});
        if (savePlots)
            saveStr = sprintf('%s/bpmV_CORR_Mon2_RefBPM2V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon2_BPMV,corrRefBPM2S_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon2 Correlation')
        ylabel('BPM CT 285S Correlation')
        title({'Correlation with BPMV: Mon2 vs. RefBPM2S', sprintf('(correl = %.2f)',bpmV_CORR_Mon2_RefBPM2S)});
        if (savePlots)
            saveStr = sprintf('%s/bpmV_CORR_Mon2_RefBPM2S',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end
    end
    
    if (useRefBPM1)
        scatter(corrMon2_BPMV,corrRefBPM1H_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon2 Correlation')
        ylabel('BPM CL 502H Correlation')
        title({'Correlation with BPMV: Mon2 vs. RefBPM1H', sprintf('(correl = %.2f)',bpmV_CORR_Mon2_RefBPM1H)});
        if (savePlots)
            saveStr = sprintf('%s/bpmV_CORR_Mon2_RefBPM1H',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon2_BPMH,corrRefBPM1V_BPMH);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon2 Correlation')
        ylabel('BPM CL 502V Correlation')
        title({'Correlation with BPMH: Mon2 vs. RefBPM1V', sprintf('(correl = %.2f)',bpmH_CORR_Mon2_RefBPM1V)});
        if (savePlots)
            saveStr = sprintf('%s/bpmH_CORR_Mon2_RefBPM1V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon2_BPMS,corrRefBPM1V_BPMS);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon2 Correlation')
        ylabel('BPM CL 502V Correlation')
        title({'Correlation with BPMS: Mon2 vs. RefBPM1V', sprintf('(correl = %.2f)',bpmS_CORR_Mon2_RefBPM1V)});
        if (savePlots)
            saveStr = sprintf('%s/bpmS_CORR_Mon2_RefBPM1V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon2_BPMV,corrRefBPM1V_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon2 Correlation')
        ylabel('BPM CL 502V Correlation')
        title({'Correlation with BPMV: Mon2 vs. RefBPM1V', sprintf('(correl = %.2f)',bpmV_CORR_Mon2_RefBPM1V)});
        if (savePlots)
            saveStr = sprintf('%s/bpmV_CORR_Mon2_RefBPM1V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon2_BPMV,corrRefBPM1S_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon2 Correlation')
        ylabel('BPM CL 502S Correlation')
        title({'Correlation with BPMV: Mon2 vs. RefBPM1S', sprintf('(correl = %.2f)',bpmV_CORR_Mon2_RefBPM1S)});
        if (savePlots)
            saveStr = sprintf('%s/bpmV_CORR_Mon2_RefBPM1S',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end
    end
end