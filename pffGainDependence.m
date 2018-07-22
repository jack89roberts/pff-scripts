rud = 0.9;
sd = 1;
su = 0.9;
g = -0.5:0.01:2.5;

rupff = sqrt((((rud^2 - 1)*sd^2)./(sd^2 + (g.^2).*su^2 - 2.*g*rud*su*sd))+1);
spff = sqrt( sd^2 + (g.^2).*(su^2) - 2.*g*rud*su*sd );

figure;
plot(g,rupff)

figure;
plot(g,spff)

%%
rud = 0.42;
sd = 2.109;
su = 1.296;
g = dataSetValues.*(rud*(sd/su))/40;

rupff = sqrt((((rud.^2 - 1).*sd.^2)./(sd.^2 + (g.^2).*su.^2 - 2.*g.*rud.*su.*sd))+1);
spff = sqrt( sd.^2 + (g.^2).*(su.^2) - 2.*g.*rud.*su.*sd );

figure;
plot(g,rupff,'o')

figure;
plot(g,spff,'o')

