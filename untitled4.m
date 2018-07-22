%%
myDir = '/home/jack/Documents/pff-thesis/Figures/feedforward';
close all;
D = dir([myDir '/*.fig']);
for i=1:length(D)
    openfig([myDir '/' D(i).name]);
    grid off;
    title('');
    matlab2tikz([myDir '/' D(i).name(1:end-3) 'tikz'],'width','0.75\textwidth');
    close;
end