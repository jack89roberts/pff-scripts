%%% Correlation between phase correlation and bpm correlation

%% Mon3

if (useCT285)
    scatter(corrMon3_BPMH,corrBPMCT285H_BPMH);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon3 Correlation')
    ylabel('BPM CT 285H Correlation')
    title({'Correlation with BPMH: Mon3 vs. CT285H', sprintf('(correl = %.2f)',bpmH_CORR_Mon3_CT285H)});
    if (savePlots)
        saveStr = sprintf('%sbpmH_CORR_Mon3_CT285H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    scatter(corrMon3_BPMS,corrBPMCT285H_BPMS);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon3 Correlation')
    ylabel('BPM CT 285H Correlation')
    title({'Correlation with BPMS: Mon3 vs. CT285H', sprintf('(correl = %.2f)',bpmS_CORR_Mon3_CT285H)});
    if (savePlots)
        saveStr = sprintf('%sbpmS_CORR_Mon3_CT285H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    scatter(corrMon3_BPMH,corrBPMCT285S_BPMH);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon3 Correlation')
    ylabel('BPM CT 285S Correlation')
    title({'Correlation with BPMH: Mon3 vs. CT285S', sprintf('(correl = %.2f)',bpmH_CORR_Mon3_CT285S)});
    if (savePlots)
        saveStr = sprintf('%sbpmH_CORR_Mon3_CT285S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    scatter(corrMon3_BPMS,corrBPMCT285S_BPMS);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon3 Correlation')
    ylabel('BPM CT 285S Correlation')
    title({'Correlation with BPMS: Mon3 vs. CT285S', sprintf('(correl = %.2f)',bpmS_CORR_Mon3_CT285S)});
    if (savePlots)
        saveStr = sprintf('%sbpmS_CORR_Mon3_CT285S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end
end

if (useCL502)
    scatter(corrMon3_BPMH,corrBPMCL502H_BPMH);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon3 Correlation')
    ylabel('BPM CL 502H Correlation')
    title({'Correlation with BPMH: Mon3 vs. CL502H', sprintf('(correl = %.2f)',bpmH_CORR_Mon3_CL502H)});
    if (savePlots)
        saveStr = sprintf('%sbpmH_CORR_Mon3_CL502H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    scatter(corrMon3_BPMS,corrBPMCL502H_BPMS);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon3 Correlation')
    ylabel('BPM CL 502H Correlation')
    title({'Correlation with BPMS: Mon3 vs. CL502H', sprintf('(correl = %.2f)',bpmS_CORR_Mon3_CL502H)});
    if (savePlots)
        saveStr = sprintf('%sbpmS_CORR_Mon3_CL502H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end


    scatter(corrMon3_BPMH,corrBPMCL502S_BPMH);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon3 Correlation')
    ylabel('BPM CL 502S Correlation')
    title({'Correlation with BPMH: Mon3 vs. CL502S', sprintf('(correl = %.2f)',bpmH_CORR_Mon3_CL502S)});
    if (savePlots)
        saveStr = sprintf('%sbpmH_CORR_Mon3_CL502S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    scatter(corrMon3_BPMS,corrBPMCL502S_BPMS);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon3 Correlation')
    ylabel('BPM CL 502S Correlation')
    title({'Correlation with BPMS: Mon3 vs. CL502S', sprintf('(correl = %.2f)',bpmS_CORR_Mon3_CL502S)});
    if (savePlots)
        saveStr = sprintf('%sbpmS_CORR_Mon3_CL502S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end
end

if (useVertical)
    if (useCT285)
        scatter(corrMon3_BPMV,corrBPMCT285H_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon3 Correlation')
        ylabel('BPM CT 285H Correlation')
        title({'Correlation with BPMV: Mon3 vs. CT285H', sprintf('(correl = %.2f)',bpmV_CORR_Mon3_CT285H)});
        if (savePlots)
            saveStr = sprintf('%sbpmV_CORR_Mon3_CT285H',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon3_BPMH,corrBPMCT285V_BPMH);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon3 Correlation')
        ylabel('BPM CT 285V Correlation')
        title({'Correlation with BPMH: Mon3 vs. CT285V', sprintf('(correl = %.2f)',bpmH_CORR_Mon3_CT285V)});
        if (savePlots)
            saveStr = sprintf('%sbpmH_CORR_Mon3_CT285V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon3_BPMS,corrBPMCT285V_BPMS);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon3 Correlation')
        ylabel('BPM CT 285V Correlation')
        title({'Correlation with BPMS: Mon3 vs. CT285V', sprintf('(correl = %.2f)',bpmS_CORR_Mon3_CT285V)});
        if (savePlots)
            saveStr = sprintf('%sbpmS_CORR_Mon3_CT285V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon3_BPMV,corrBPMCT285V_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon3 Correlation')
        ylabel('BPM CT 285V Correlation')
        title({'Correlation with BPMV: Mon3 vs. CT285V', sprintf('(correl = %.2f)',bpmV_CORR_Mon3_CT285V)});
        if (savePlots)
            saveStr = sprintf('%sbpmV_CORR_Mon3_CT285V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon3_BPMV,corrBPMCT285S_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon3 Correlation')
        ylabel('BPM CT 285S Correlation')
        title({'Correlation with BPMV: Mon3 vs. CT285S', sprintf('(correl = %.2f)',bpmV_CORR_Mon3_CT285S)});
        if (savePlots)
            saveStr = sprintf('%sbpmV_CORR_Mon3_CT285S',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end
    end
    
    if (useCL502)
        scatter(corrMon3_BPMV,corrBPMCL502H_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon3 Correlation')
        ylabel('BPM CL 502H Correlation')
        title({'Correlation with BPMV: Mon3 vs. CL502H', sprintf('(correl = %.2f)',bpmV_CORR_Mon3_CL502H)});
        if (savePlots)
            saveStr = sprintf('%sbpmV_CORR_Mon3_CL502H',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon3_BPMH,corrBPMCL502V_BPMH);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon3 Correlation')
        ylabel('BPM CL 502V Correlation')
        title({'Correlation with BPMH: Mon3 vs. CL502V', sprintf('(correl = %.2f)',bpmH_CORR_Mon3_CL502V)});
        if (savePlots)
            saveStr = sprintf('%sbpmH_CORR_Mon3_CL502V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon3_BPMS,corrBPMCL502V_BPMS);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon3 Correlation')
        ylabel('BPM CL 502V Correlation')
        title({'Correlation with BPMS: Mon3 vs. CL502V', sprintf('(correl = %.2f)',bpmS_CORR_Mon3_CL502V)});
        if (savePlots)
            saveStr = sprintf('%sbpmS_CORR_Mon3_CL502V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon3_BPMV,corrBPMCL502V_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon3 Correlation')
        ylabel('BPM CL 502V Correlation')
        title({'Correlation with BPMV: Mon3 vs. CL502V', sprintf('(correl = %.2f)',bpmV_CORR_Mon3_CL502V)});
        if (savePlots)
            saveStr = sprintf('%sbpmV_CORR_Mon3_CL502V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon3_BPMV,corrBPMCL502S_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon3 Correlation')
        ylabel('BPM CL 502S Correlation')
        title({'Correlation with BPMV: Mon3 vs. CL502S', sprintf('(correl = %.2f)',bpmV_CORR_Mon3_CL502S)});
        if (savePlots)
            saveStr = sprintf('%sbpmV_CORR_Mon3_CL502S',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end
    end

end

%% Mon2

if (useCT285)
    scatter(corrMon2_BPMH,corrBPMCT285H_BPMH);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon2 Correlation')
    ylabel('BPM CT 285H Correlation')
    title({'Correlation with BPMH: Mon2 vs. CT285H', sprintf('(correl = %.2f)',bpmH_CORR_Mon2_CT285H)});
    if (savePlots)
        saveStr = sprintf('%sbpmH_CORR_Mon2_CT285H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    scatter(corrMon2_BPMS,corrBPMCT285H_BPMS);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon2 Correlation')
    ylabel('BPM CT 285H Correlation')
    title({'Correlation with BPMS: Mon2 vs. CT285H', sprintf('(correl = %.2f)',bpmS_CORR_Mon2_CT285H)});
    if (savePlots)
        saveStr = sprintf('%sbpmS_CORR_Mon2_CT285H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    scatter(corrMon2_BPMH,corrBPMCT285S_BPMH);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon2 Correlation')
    ylabel('BPM CT 285S Correlation')
    title({'Correlation with BPMH: Mon2 vs. CT285S', sprintf('(correl = %.2f)',bpmH_CORR_Mon2_CT285S)});
    if (savePlots)
        saveStr = sprintf('%sbpmH_CORR_Mon2_CT285S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    scatter(corrMon2_BPMS,corrBPMCT285S_BPMS);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon2 Correlation')
    ylabel('BPM CT 285S Correlation')
    title({'Correlation with BPMS: Mon2 vs. CT285S', sprintf('(correl = %.2f)',bpmS_CORR_Mon2_CT285S)});
    if (savePlots)
        saveStr = sprintf('%sbpmS_CORR_Mon2_CT285S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end
end

if (useCL502)
    scatter(corrMon2_BPMH,corrBPMCL502H_BPMH);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon2 Correlation')
    ylabel('BPM CL 502H Correlation')
    title({'Correlation with BPMH: Mon2 vs. CL502H', sprintf('(correl = %.2f)',bpmH_CORR_Mon2_CL502H)});
    if (savePlots)
        saveStr = sprintf('%sbpmH_CORR_Mon2_CL502H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    scatter(corrMon2_BPMS,corrBPMCL502H_BPMS);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon2 Correlation')
    ylabel('BPM CL 502H Correlation')
    title({'Correlation with BPMS: Mon2 vs. CL502H', sprintf('(correl = %.2f)',bpmS_CORR_Mon2_CL502H)});
    if (savePlots)
        saveStr = sprintf('%sbpmS_CORR_Mon2_CL502H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end


    scatter(corrMon2_BPMH,corrBPMCL502S_BPMH);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon2 Correlation')
    ylabel('BPM CL 502S Correlation')
    title({'Correlation with BPMH: Mon2 vs. CL502S', sprintf('(correl = %.2f)',bpmH_CORR_Mon2_CL502S)});
    if (savePlots)
        saveStr = sprintf('%sbpmH_CORR_Mon2_CL502S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    scatter(corrMon2_BPMS,corrBPMCL502S_BPMS);
    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)
    grid on;
    xlabel('Mon2 Correlation')
    ylabel('BPM CL 502S Correlation')
    title({'Correlation with BPMS: Mon2 vs. CL502S', sprintf('(correl = %.2f)',bpmS_CORR_Mon2_CL502S)});
    if (savePlots)
        saveStr = sprintf('%sbpmS_CORR_Mon2_CL502S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end
end

if (useVertical)
    if (useCT285)
        scatter(corrMon2_BPMV,corrBPMCT285H_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon2 Correlation')
        ylabel('BPM CT 285H Correlation')
        title({'Correlation with BPMV: Mon2 vs. CT285H', sprintf('(correl = %.2f)',bpmV_CORR_Mon2_CT285H)});
        if (savePlots)
            saveStr = sprintf('%sbpmV_CORR_Mon2_CT285H',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon2_BPMH,corrBPMCT285V_BPMH);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon2 Correlation')
        ylabel('BPM CT 285V Correlation')
        title({'Correlation with BPMH: Mon2 vs. CT285V', sprintf('(correl = %.2f)',bpmH_CORR_Mon2_CT285V)});
        if (savePlots)
            saveStr = sprintf('%sbpmH_CORR_Mon2_CT285V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon2_BPMS,corrBPMCT285V_BPMS);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon2 Correlation')
        ylabel('BPM CT 285V Correlation')
        title({'Correlation with BPMS: Mon2 vs. CT285V', sprintf('(correl = %.2f)',bpmS_CORR_Mon2_CT285V)});
        if (savePlots)
            saveStr = sprintf('%sbpmS_CORR_Mon2_CT285V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon2_BPMV,corrBPMCT285V_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon2 Correlation')
        ylabel('BPM CT 285V Correlation')
        title({'Correlation with BPMV: Mon2 vs. CT285V', sprintf('(correl = %.2f)',bpmV_CORR_Mon2_CT285V)});
        if (savePlots)
            saveStr = sprintf('%sbpmV_CORR_Mon2_CT285V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon2_BPMV,corrBPMCT285S_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon2 Correlation')
        ylabel('BPM CT 285S Correlation')
        title({'Correlation with BPMV: Mon2 vs. CT285S', sprintf('(correl = %.2f)',bpmV_CORR_Mon2_CT285S)});
        if (savePlots)
            saveStr = sprintf('%sbpmV_CORR_Mon2_CT285S',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end
    end
    
    if (useCL502)
        scatter(corrMon2_BPMV,corrBPMCL502H_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon2 Correlation')
        ylabel('BPM CL 502H Correlation')
        title({'Correlation with BPMV: Mon2 vs. CL502H', sprintf('(correl = %.2f)',bpmV_CORR_Mon2_CL502H)});
        if (savePlots)
            saveStr = sprintf('%sbpmV_CORR_Mon2_CL502H',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon2_BPMH,corrBPMCL502V_BPMH);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon2 Correlation')
        ylabel('BPM CL 502V Correlation')
        title({'Correlation with BPMH: Mon2 vs. CL502V', sprintf('(correl = %.2f)',bpmH_CORR_Mon2_CL502V)});
        if (savePlots)
            saveStr = sprintf('%sbpmH_CORR_Mon2_CL502V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon2_BPMS,corrBPMCL502V_BPMS);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon2 Correlation')
        ylabel('BPM CL 502V Correlation')
        title({'Correlation with BPMS: Mon2 vs. CL502V', sprintf('(correl = %.2f)',bpmS_CORR_Mon2_CL502V)});
        if (savePlots)
            saveStr = sprintf('%sbpmS_CORR_Mon2_CL502V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon2_BPMV,corrBPMCL502V_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon2 Correlation')
        ylabel('BPM CL 502V Correlation')
        title({'Correlation with BPMV: Mon2 vs. CL502V', sprintf('(correl = %.2f)',bpmV_CORR_Mon2_CL502V)});
        if (savePlots)
            saveStr = sprintf('%sbpmV_CORR_Mon2_CL502V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        scatter(corrMon2_BPMV,corrBPMCL502S_BPMV);
        set(gca,'FontSize',16)
        set(findall(gcf,'type','text'),'FontSize',16)
        grid on;
        xlabel('Mon2 Correlation')
        ylabel('BPM CL 502S Correlation')
        title({'Correlation with BPMV: Mon2 vs. CL502S', sprintf('(correl = %.2f)',bpmV_CORR_Mon2_CL502S)});
        if (savePlots)
            saveStr = sprintf('%sbpmV_CORR_Mon2_CL502S',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end
    end
end