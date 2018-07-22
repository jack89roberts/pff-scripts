%function phase = getPhaseMixerDiode( Mixer, Diode, MixerSensitivity, DiodeSensitivity, CalibrationConstant )
function phase = getPhaseMixerDiode( Mixer, Diode, CalibrationAmplitude, CalibrationOffset, useDiode )
% get phase (degrees) from Mixer and Diode channel inputs
%   for frascati monitors
%   CalibrationOffset: zero if not given.
%   useDiode: 0=justMixer, 1=mixer/sqrt(Diode), 2=mixer/Diode. Default: 0
%   if diode empty, 1 otherwise.
    
    if nargin<5
        if (isempty(Diode))
            useDiode = 0;
        else
            useDiode = 1;
        end
    end
    
    if nargin<4
        CalibrationOffset = 0;
    end

    if (useDiode==1) % Mixer/sqrt(Diode)
        phase = real(asind( ((Mixer ./ sqrt(abs(Diode))) - CalibrationOffset) ./CalibrationAmplitude ));
    
    elseif (useDiode==2) % Mixer/Diode
        phase = real(asind( ((Mixer ./ abs(Diode)) - CalibrationOffset) ./CalibrationAmplitude ));    
    
    else % Mixer only
        phase = real(asind( (Mixer - CalibrationOffset) ./CalibrationAmplitude ));        
    end
    
    
end

