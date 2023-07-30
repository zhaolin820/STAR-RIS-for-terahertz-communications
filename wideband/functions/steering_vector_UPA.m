function [a] = steering_vector_UPA(theta_h, theta_v, N_h, N_v, fm, fc)
%Calculate steering vector for UPAs
%  [a] = steering_vector_UPA(theta_h, theta_v, N_h, N_v)
%Inputs:
%   theta_h: azimuth angle
%   theta_v: elevation angle
%   N_h x N_v: number of antennas
%   f_m: carrier frequency
%   f_c: central frequency
%Outputs:
%   a: steering vector
%Date: 27/07/2023
%Author: Zhaolin Wang

h = 0:(N_h-1);
a_h = exp(1i * pi * h * fm/fc * sin(theta_h) * sin(theta_v)).';

v = 0:(N_v-1);
a_v = exp(1i * pi * v * fm/fc * cos(theta_v)).';

a = kron(a_h, a_v);
end


