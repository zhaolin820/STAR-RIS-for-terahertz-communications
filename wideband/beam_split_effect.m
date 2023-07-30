clc
clear all
close all
addpath("functions\");
para = para_init();

theta = 45; % target direction
theta_all = theta-8:0.01:theta+8;
theta_l = sin(theta*pi/180);


%% TTD-based hybrid beamforming
N_sub = para.N / para.N_T;
f = steering_vector_ULA(theta*pi/180, para.N, para.fc, para.fc);
F_RF = zeros(para.N, para.N_T);
t = zeros(para.N_T, 1);

for i = 1:para.N_T
    F_RF((i-1)*N_sub+1:i*N_sub, i) = f((i-1)*N_sub+1:i*N_sub) * exp(-1i*pi*(i-1)*N_sub*theta_l);

    if theta_l >= 0
        t(i) = (i-1) * (N_sub*theta_l/2) * 1/para.fc;
    else
        t(i) = (para.N_T-1) * abs(N_sub*theta_l/2) * 1/para.fc + (i-1)* (N_sub*theta_l/2) * 1/para.fc;
    end
end

figure; 
subplot(2,1,2); hold on; box on;

% central frequency
fm = F_RF*exp(1i*2*pi*para.fc*t);
array_gain_c = zeros(length(theta_all), 1);
for i = 1:length(theta_all)
    a = steering_vector_ULA(theta_all(i)*pi/180, para.N, para.fc, para.fc);
    array_gain_c(i) = abs(fm' * a);
end
array_gain_c = array_gain_c / max(array_gain_c);
plot(theta_all, array_gain_c, 'b');

% maximum frequency
fm = F_RF*exp(1i*2*pi*para.fm_all(end)*t);
array_gain_Mc = zeros(length(theta_all), 1);
for i = 1:length(theta_all)
    a = steering_vector_ULA(theta_all(i)*pi/180, para.N, para.fm_all(end), para.fc);
    array_gain_Mc(i) = abs(fm' * a);
end
array_gain_Mc = array_gain_Mc / max(array_gain_Mc);
plot(theta_all, array_gain_Mc, 'r');

% minimum frequency
fm = F_RF*exp(1i*2*pi*para.fm_all(1)*t);
array_gain_1 = zeros(length(theta_all), 1);
for i = 1:length(theta_all)
    a = steering_vector_ULA(theta_all(i)*pi/180, para.N, para.fm_all(1), para.fc);
    array_gain_1(i) = abs(fm' * a);
end
array_gain_1 = array_gain_1 / max(array_gain_1);
plot(theta_all, array_gain_1, 'g');

xlabel("Angle (degree)");
ylabel("Normalized array gain");
legend("$f_c$", "$f_1$", "$f_{M_c}$", 'Interpreter','latex');


%% Conventional hybrid beamforming

f = steering_vector_ULA(theta*pi/180, para.N, para.fc, para.fc);

subplot(2,1,1); hold on; box on;

% central frequency
array_gain_c = zeros(length(theta_all), 1);
for i = 1:length(theta_all)
    a = steering_vector_ULA(theta_all(i)*pi/180, para.N, para.fc, para.fc);
    array_gain_c(i) = abs(f' * a);
end
array_gain_c = array_gain_c / max(array_gain_c);
plot(theta_all, array_gain_c, 'b');

% maximum frequency
array_gain_Mc = zeros(length(theta_all), 1);
for i = 1:length(theta_all)
    a = steering_vector_ULA(theta_all(i)*pi/180, para.N, para.fm_all(end), para.fc);
    array_gain_Mc(i) = abs(f' * a);
end
array_gain_Mc = array_gain_Mc / max(array_gain_Mc);
plot(theta_all, array_gain_Mc, 'r');

% minimum frequency
array_gain_1 = zeros(length(theta_all), 1);
for i = 1:length(theta_all)
    a = steering_vector_ULA(theta_all(i)*pi/180, para.N, para.fm_all(1), para.fc);
    array_gain_1(i) = abs(f' * a);
end
array_gain_1 = array_gain_1 / max(array_gain_1);
plot(theta_all, array_gain_1, 'g');

xlabel("Angle (degree)");
ylabel("Normalized array gain");
legend("$f_c$", "$f_1$", "$f_{M_c}$", 'Interpreter','latex');