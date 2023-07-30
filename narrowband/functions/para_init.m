function [para] = para_init()
%Construct a struct of the initial values for all the parameters 
%  [values] = para_init()
%Inputs:
%   None
%Outputs:
%   values: a struct
%Date: 27/06/2021
%Author: Zhaolin Wang

para.f = 1e11; % carrier frequency 0.1 THz
B = 1e8; % 100 MHz system bandwidth 

para.noise_dB = -174; % noise power in dBm/Hz
para.noise_dB = para.noise_dB + 10*log10(B) - 30; % noise power in dBW

para.Gt = 25; % dBi, transmit antenna gain
para.Gr = 15; % dBi, receive antenna gain

para.N = 128; % number of antennas
para.N_RF = 4; % number of RF chains

para.M_h = 6;
para.M_v = 6;
para.M = para.M_h*para.M_v; % reflecting elements at RIS

para.Pt = 10^(20/10); % W overall transmit power

para.L = 4; % number of paths between BS and STARS
para.Lk = 4; % number of paths between STARS and each user
para.K = 4; % user number

P_UE = 100; % 100 mW, user power 
P_STAR_independent = para.M * 1/2*(7+2*8)*0.33 + 1e4; % mW, STAR-RIS power
P_STAR_coupled = para.M * 1/2*(7+8+1)*0.33+ 1e4; % mW, STAR-RIS power
P_RIS = para.M * 1/2*(8)*0.33+ 1e4;
P_BB = 300; % 300 mW, baseband processing power 
P_RF = 200; % 200 mW, RF chain power
P_PS = 10; % 20 mW, PS power
P_BS = 3e3; % 3 W, BS power
para.Pc_HB_idp =1e-3 * (P_BS + P_BB + para.N_RF * P_RF + para.N_RF * para.N * P_PS + P_STAR_independent + para.K*P_UE); 
para.Pc_HB_coup =1e-3 * (P_BS + P_BB + para.N_RF * P_RF + para.N_RF * para.N * P_PS + P_STAR_coupled + para.K*P_UE); 

para.xi = 1e-3 * 100; % 100 mW/(bit/s/Hz), rate-dependent power consumption factor.


end

