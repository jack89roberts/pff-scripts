%% improved zero crossing calculating and error bar

myZeroes = -calibrationFactors(:,3)./calibrationFactors(:,2);
n = zeros(1,nMons);
for mon=1:nMons
    dydxAtZero = calibrationFactors(mon,1)*calibrationFactors(mon,2)*cos(calibrationFactors(mon,2)*myZeroes(mon) + calibrationFactors(mon,3)); % Ab*cos(bx+c)

    while (~isnan(myZeroes(mon)) && (myZeroes(mon)<0 || dydxAtZero>0))
        myZeroes(mon) = myZeroes(mon) + pi/calibrationFactors(mon,2);
        n(mon) = n(mon)+1;
        dydxAtZero = calibrationFactors(mon,1)*calibrationFactors(mon,2)*cos(calibrationFactors(mon,2)*myZeroes(mon) + calibrationFactors(mon,3)); % Ab*cos(bx+c)
    end
end

myZeroes_err = NaN(1,nMons);
for mon=1:nMons
    myZeroes_err(mon) = (calConstErr(mon,3)./calibrationFactors(mon,2)).^2;
    myZeroes_err(mon) = myZeroes_err(mon) + ( (calConstErr(mon,2)./(calibrationFactors(mon,2).^2)).*(calibrationFactors(mon,3)-n(mon)*pi) ).^2;
    myZeroes_err(mon) = sqrt(myZeroes_err(mon));
end

myZeroes'
myZeroes_err

%%
% [x,y,n,d]=getDataFromFig
m3 = y{2};
m2 = y{4};
m1 = y{6};
x = x{2};
A = [1.1167 1.064 0.99];
b = [0.0508 0.0474 0.0513];
c = [2.29 1.36 4.02];
d = [0.086 0.069 0.15];

fitM1 = A(1).*sin(b(1).*x + c(1)) + d(1);
fitM2 = A(2).*sin(b(2).*x + c(2)) + d(2);
fitM3 = A(3).*sin(b(3).*x + c(3)) + d(3);

figure;
plot(x,asind((m1-fitM1)./A(1)),'b');
hold all
plot(x,asind((m2-fitM2)./A(2)),'r');
plot(x,asind((m3-fitM3)./A(3)),'g');
