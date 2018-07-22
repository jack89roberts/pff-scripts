%%
myDir = '/home/jack/Documents/pff-thesis/Figures/feedforward';
close all;
D = dir([myDir '/*.fig']);
for i=1:length(D)
    fprintf('%d of %d: %s... ', i,length(D),D(i).name);
    try
        openfig([myDir '/' D(i).name]);
        grid off;
        title('');
        legend off;
        matlab2tikz([myDir '/' D(i).name(1:end-3) 'tikz'],'width','0.75\textwidth','height','0.555\textwidth','showInfo',false);
        fprintf ('done!\n');
    catch
        fprintf('failed to convert!\n');
    end
    close;
end