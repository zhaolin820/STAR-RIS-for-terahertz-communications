function [loss] = getAbsLoss(f, d, HITRANparams)
%Calculate the absorption loss
%  [loss] = getAbsLoss(f, d, HITRANparams)
%Inputs:
%   f: frequency
%   d: distance
%   HITRANparams: HITRAN database
%Outputs:
%   loss: absorption loss
%Author: Michele Polese
%Link: https://github.com/mychele/toward-e2e-6g-terahertz-networks

[minDiff, closestFreqIndex] = min(abs(HITRANparams(:, 1) - f));

kfParam = 0;
if (abs(minDiff) < 9.894e8) % this 9.894e8 is from the 
	kfParam = HITRANparams(closestFreqIndex, 2);
end

loss = kfParam * d * 10 * log10(exp(1));

end

