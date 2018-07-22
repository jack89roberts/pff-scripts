function [ mon3Phases ] = subtractChanOffset( mon3Phases, chanOffsets, mon1A, mon1D, refOffset,polarity )
%subtractChanOffset Takes downstream phases and subtracts channel offset
%from data to take in to account downstream phase shifts resulting from
%differences in PFF setup during a dataset.
%   [ mon3Phases ] = subtractChanOffset( mon3Phases, chanOffsets, mon1A, mon1D, refOffset, polarity )
%
%   mon3Phases - mon3 sample by sample phases (pulses,samples)
%
%   chanOffsets - channel offset at each pulse (from FONT data)
%
%   mon1A, mon1D - Mon1 calibration constants, to convert channel offset in
%   to expected downstream phase shift.
%
%   refOffset - which phase offset to correct to. (default: mean of
%   phaseOffsets).
%
%   polarity - whether to add or subtract (default: 1).

    if (nargin < 6 || isempty(polarity))
        polarity = 1;
    end
    
    phaseOffsets = polarity.*asind((chanOffsets-mon1D)./mon1A);
    if (nargin < 5 || isempty(refOffset))
        phaseOffsets = phaseOffsets-nanmean(phaseOffsets);
    else
        phaseOffsets = phaseOffsets-refOffset;
    end

    nPulses = length(chanOffsets);
    for p=1:nPulses
        mon3Phases(p,:) = mon3Phases(p,:) + phaseOffsets(p);
    end

end

