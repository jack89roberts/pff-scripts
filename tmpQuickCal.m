figure;
mix2 = squeeze(mixers(3,:,200));
plot(mix2');

x = 1:length(mix2);
y = mix2;
x = x(1:100);
y = y(1:100);

phase = getPhaseMixerDiode(squeeze(mixers(2,:,:)),[],0.325,0.025);