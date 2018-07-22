%%
myDir = '/home/jack/Documents/pff-thesis/Figures/propagation';
close all;
D = dir([myDir '/*.fig']);
for i=1:length(D)
    fprintf('%d of %d: %s... ', i,length(D),D(i).name);
    try
        openfig([myDir '/' D(i).name]);
        grid off;
        title('');
        legend off;
        pos = get(gcf,'Position');
        pos(4) = 0.9*pos(4);
        set(gcf,'Position',pos);
        %matlab2tikz([myDir '/' D(i).name(1:end-3) 'tikz'],'width','0.75\textwidth','height','0.555\textwidth','showInfo',false);
        savePlot(myDir, D(i).name(1:end-4),[0,1,0,1]); % saves fig and pdf
        fprintf ('done!\n');
    catch
        fprintf('failed to convert!\n');
    end
    close;
end