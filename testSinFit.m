x = linspace(1,360);
y = sind(x);
plot(x,y,'o');

dydx = NaN*zeros(1,length(x)-1);
for i=1:(length(x)-1)
    dydx(i) = (y(i+1)-y(i))./(x(i+1)-x(i));
end

%plot(x(2:length(x)),dydx)

plot(x(34:67),y(34:67));

fit(x,y,'poly1')