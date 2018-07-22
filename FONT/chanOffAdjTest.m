% close all;

N = 1000;
avg = 100;

beam = (1:N) + random('norm',0,50,1,N);
beam = [beam (1000*ones(1,1000)+ random('norm',0,50,1,N)) ((N:-1:1)+ random('norm',0,50,1,N)) (zeros(1,1000)+ random('norm',0,50,1,N))];
% beam = random('norm',0,50,1,N) + 1500;

N = length(beam);

off = NaN(1,N);
adc = NaN(1,N);

off(1) = 0;
adc(1) = beam(1);

for i=2:N
    adc(i) = beam(i) + off(i-1);
    
    % running
    if i>avg+1
        off(i) = nanmean(off((i-1-avg):(i-1))) - nanmean(adc((i-avg):i));
    else
        off(i) = nanmean(off(1:(i-1))) - nanmean(adc(1:i));
    end
    
    % box car
%     if (mod(i,avg)==0)
%         off(i) = off(i-1) - nanmean(adc((i-avg+1):i));
%     else
%         off(i) = off(i-1);
%     end
end

figure;
plot(beam);
title('beam');

figure;
plot(adc)
title('adc')

figure;
plot(off)
title('off')
