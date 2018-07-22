%% Input parameters

relativeSign = 1;
testMode = false;

%
% wiggle mode: change phase in steps of stepSize between +/- maxStep.
%
setting = 'wiggle';
stepSize = 0.1;
maxStep = 0.3;
stepSign = 1;

%
% step mode: change phase by +/- stepSize.
%
% setting = 'step';
% stepSize = 0.5;
% stepSign = 1;

%
% gauss mode: change phase with gaussian distribution, std = stepStd, don't
% allow phase shifts above maxStep.
%
% setting = 'gauss';
% stepStd = 0.3;
% maxStep = 1.0;

%% wiggle phase

selector = 'SCT.USER.SETUP';
%Setting/PhaseSh_SP
initPhase2 = matlabJapc.staticGetSignal(selector, 'CK.LL-MKS02/Acquisition#PhaseSH_VAL');
initPhase3 = matlabJapc.staticGetSignal(selector, 'CK.LL-MKS03/Acquisition#PhaseSH_VAL');
currPhase2 = initPhase2;
currPhase3 = initPhase3;

lengthHistory = 100;
phaseHistory2 = NaN*ones(1,lengthHistory);
phaseHistory3 = NaN*ones(1,lengthHistory);
pulseIndex = 1;

if (testMode)
    fprintf(1,'----------------------------------------------------\n');
    fprintf(1,'TEST MODE ENABLED - NEW PHASE VALUES WILL NOT BE SET\n');
    fprintf(1,'----------------------------------------------------\n');
end

fprintf(1, 'Starting MKS02 and MKS03 phase shifts in %s mode...\n',setting);
fprintf(1, 'Ctrl-C to stop, then remember to set init phases afterwards (last section in script)\n');

while (true)
    if (strcmp(setting,'wiggle'))
        if abs(currPhase2 - initPhase2) >= maxStep
            stepSign = -stepSign;
        end
        currPhase2 = currPhase2 + stepSign.*stepSize;
        currPhase3 = currPhase3 + relativeSign.*stepSign.*stepSize;
    
    elseif (strcmp(setting,'step'))
        if (pulseIndex == 1)
            currPhase2 = currPhase2 + stepSign.*stepSize;
            currPhase3 = currPhase3 + relativeSign.*stepSign.*stepSize;
        else
            currPhase2 = currPhase2 + 2.*stepSign.*stepSize;
            currPhase3 = currPhase3 + 2.*relativeSign.*stepSign.*stepSize; 
        end
        stepSign = -stepSign;
    elseif (strcmp(setting,'gauss'))
        stepSize = maxStep + 1;
        while (abs(stepSize) > maxStep)
            stepSize = random('norm',0,stepStd);
        end
        currPhase2 = initPhase2 + stepSize;
        currPhase3 = initPhase3 + relativeSign.*stepSize;
    end
    
    if (~testMode)
        matlabJapc.staticSetSignal(selector, 'CK.LL-MKS02/Setting#PhaseSh_SP', currPhase2);
        matlabJapc.staticSetSignal(selector, 'CK.LL-MKS03/Setting#PhaseSh_SP', currPhase3);
    end
    
    
    if (pulseIndex > lengthHistory)
        phaseHistory2 =  circshift(phaseHistory2,[0 -1]);
        phaseHistory3 =  circshift(phaseHistory3,[0 -1]);
        phaseHistory2(lengthHistory) = currPhase2; 
        phaseHistory3(lengthHistory) = currPhase3;
        pulseIndex = pulseIndex + 1;
    else
        phaseHistory2(pulseIndex) = currPhase2; 
        phaseHistory3(pulseIndex) = currPhase3;
        pulseIndex = pulseIndex + 1;
    end
    plot(1:lengthHistory,phaseHistory2-initPhase2,'bo-');
    hold all;
    plot(1:lengthHistory,phaseHistory3-initPhase3,'ro-');
    if (testMode)
        title('Simulated Values To Be Sent (TEST MODE)');
    else
        title('Values Sent to MKS02 and MKS03');
    end
    xlabel('Pulse Index');
    ylabel('PhShifter Setting');
    hold off;
    
    pause(1/1.2);
end

%% reset init phases   
matlabJapc.staticSetSignal(selector, 'CK.LL-MKS02/Setting#PhaseSh_SP', initPhase2);
matlabJapc.staticSetSignal(selector, 'CK.LL-MKS03/Setting#PhaseSh_SP', initPhase3);


