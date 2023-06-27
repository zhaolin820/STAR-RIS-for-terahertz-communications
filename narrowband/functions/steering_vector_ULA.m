function [a] = steering_vector_ULA(theta, N)
%Calculate steering vector for ULAs
%  [a] = steering_vector_ULA(theta, N)
%Inputs:
%   theta: angle
%   N: number of antennas
%Outputs:
%   a: steering vector
%Date: 27/06/2023
%Author: Zhaolin Wang

n = 0:(N-1);
a = exp(1i * pi * n * sin(theta)).';
end

