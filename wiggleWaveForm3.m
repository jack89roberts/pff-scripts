close all;
%% Inputs

testMode = true;
timeShift = true;

startWiggleTime = 5500;

% shifting phase of wiggle
% mode = 'wiggle';
% stepSize = 0.4;
% maxStep = 1.6;

% shifting amplitude of wiggle
mode = 'wiggleStep';
amplitudeStepSize = 0.1;
wiggleStepSize = 0.4;
maxStep = 1.0;

%mode = 'step';
% stepSize = 0.5;
% stepFreq = 5;

minPhase = -180;
maxPhase = 180;


%% Acquire waveform
addpath('./..');
addpath('./../Klystrons');
[ctfmodroot, dataout, madxp] = ctfsettings();
addpath([ctfmodroot '/MatLab/JavaCoInterface']);
addpath([ctfmodroot '/MatLab/common']);


ampler='CK.SVPKI03A';
rpsi  ='CK.RPSI03A';
wfgenerator='CK.WFGMKS03';

profsize = JGetCoValueFesa(wfgenerator,'Profile','profileSize');
wform    = JGetCoArrayFesa(wfgenerator,'Profile','profile');
wform_orig    = wform;

for i=1:profsize
 
 t(i,1)=wform(2*i-1);
 e(i,1)=wform(2*i  );

fprintf(1,'%d t=%f e=%f \n',i,t(i,1),e(i,1));
 
  
end;
%% Calculate new waveform

if (timeShift)
    disp('Setting waveform with time shift...');
    disp('Ctrl-C to stop, remember to set initial waveform back afterwards...');
    if(testMode)
        disp('Test mode, new wave forms will not be set.');
    end
    
    figure(2)
    startPhase = 0;
    startStepSign = 1;
    startOutsideRange = false;
    firstRun = true;
    startOffset = 0;

    while (true)
        newE = e(t>startWiggleTime);
        
        if(strcmp(mode,'wiggleStep'))
            if (firstRun)
                stepSign = 1.0;
                currPhase = startPhase;
                phasesToAdd = zeros(length(newE),1);
                for i=1:length(newE)
                    phasesToAdd(i) = phasesToAdd(i) + currPhase;

                    if abs(currPhase) >= maxStep
                        stepSign = -stepSign;
                    end
                    currPhase = currPhase + stepSign.*wiggleStepSize;

                end
                startPhase = phasesToAdd./maxStep;
                wiggleAmplitude = maxStep;
                wiggleStepSign = -1;
                firstRun = false;
            else
                wiggleAmplitude = wiggleAmplitude + wiggleStepSign.*amplitudeStepSize;
                if (abs(wiggleAmplitude) >= maxStep)
                    wiggleStepSign = -wiggleStepSign;
                end
                phasesToAdd = wiggleAmplitude.*startPhase; 
            end
            newE = newE + phasesToAdd;
            newE(newE<minPhase) = minPhase;
            newE(newE>maxPhase) = maxPhase;
            
            newE = [e(t<=startWiggleTime); newE];
        
        elseif (strcmp(mode,'wiggle'))
            
            if (firstRun)
                stepSign = 1.0;
                currPhase = startPhase;
                phasesToAdd = zeros(length(newE),1);
                lastOutsideRange = false;
                for i=1:length(newE)
                    phasesToAdd(i) = phasesToAdd(i) + currPhase;

                    if abs(currPhase) >= maxStep
                        stepSign = -stepSign;
                    end
                    currPhase = currPhase + stepSign.*stepSize;

                end
                firstRun = false;
            else
               phasesToAdd = circshift(phasesToAdd,[1 0]); 
            end
            newE = newE + phasesToAdd;

            if abs(startPhase) >= maxStep
                if (~startOutsideRange)
                    startStepSign = -startStepSign;
                end
                startOutsideRange = true;
            else
                startOutsideRange = false;
            end
            startPhase = startPhase + startStepSign.*stepSize;
                
            newE(newE<minPhase) = minPhase;
            newE(newE>maxPhase) = maxPhase;
            newE = [e(t<=startWiggleTime); newE];
            
        elseif (strcmp(mode,'step'))
            
            stepSign = startStepSign;
            phasesToAdd = zeros(length(newE),1);
            for i=1:length(newE)
                phasesToAdd(i) = phasesToAdd(i) + stepSign.*stepSize;
                
                if (mod(i+startOffset,stepFreq)==0)
                    stepSign = -stepSign;
                end
            end
            
            newE = newE + phasesToAdd;
            
            newE(newE<minPhase) = minPhase;
            newE(newE>maxPhase) = maxPhase;
            newE = [e(t<=startWiggleTime); newE];
            
            if (startOffset == stepFreq-1);
                startStepSign = -startStepSign;
                startOffset = 0;
            else
                startOffset = startOffset + 1;
            end
        end
        %startOffset
        subplot(3,2,1:4)
        plot(t(t>startWiggleTime),e(t>startWiggleTime));
        hold all;
        plot(t(t>startWiggleTime),newE(t>startWiggleTime)); 
        xlim([startWiggleTime 6600])
        hold off;
        legend('Old Waveform', 'New Waveform','Location','best');
        subplot(3,2,5:6)
        plot(t(t>startWiggleTime),newE(t>startWiggleTime)-e(t>startWiggleTime))
        title('Difference Waveform');
        
        try 
            ylim([-1.1*maxStep, 1.1*maxStep]);
        catch
            ylim([-1.1*stepSize,1.1*stepSize]);
        end
        xlim([startWiggleTime 6600])
       
        
        z = wform;
        for i=1:profsize
          z(2*i - 1) = t(i,1);
          z(2*i    ) = newE(i,1);
        end

        if (~testMode)        
            JSetCoArrayFesa(wfgenerator,'Profile','profile',z);

            disp('SETTING NEW WAVE FORM...');
            fprintf(1,'Sending profile\n');
            JSetCoArrayFesa(wfgenerator,'Profile','profile',z);
            
        else
           
        end

        pause(1/1.2);
    end
    
    
else
    disp('Setting new waveform once...');
    newE = e(t>startWiggleTime);

    if (strcmp(mode,'wiggle'))
        currPhase = 0;
        stepSign = 1.0;

        for i=1:length(newE)
            newE(i) = newE(i) + currPhase;

            if abs(currPhase) >= maxStep
                stepSign = -stepSign;
            end
            currPhase = currPhase + stepSign.*stepSize;

        end

        newE = [e(t<=startWiggleTime); newE];
    end
    
    
    figure(2)
    subplot(3,2,1:4)
    plot(t,e);
    hold all;
    plot(t,newE);
    legend('Old Waveform', 'New Waveform')
    hold off;
    subplot(3,2,5:6)
    plot(t,newE-e)
    title('Difference Waveform');

    z = wform;
    for i=1:profsize
      %z(2*i - 1) = t(i,1);
      z(2*i    ) = newE(i,1);

    end
    
    
    if (~testMode)
        disp('SETTING NEW WAVE FORM...');
        fprintf(1,'Sending profile\n');
        JSetCoArrayFesa(wfgenerator,'Profile','profile',z);
    else
        disp('Test mode, wave form not set.');
    end
    

end

%% Set back initial waveform
i = input('press any key to resend the old one');
JSetCoArrayFesa(wfgenerator,'Profile','profile',wform_orig);