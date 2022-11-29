%% Channel Equalization using Least Mean Square (LMS) algorithm
% Author: SHUJAAT KHAN
clc;
clear all;
close all;
%% Channel and noise level
h = [0.9 0.3 0.5 -0.1]; % Channel
SNRr = 30;              % Noise Level
%% LMS parameters
runs = 100;        % Number of independent runs (training repetation)
eta = 5e-3;         % Learning rate / step size
order=12;           % Order of the equalizer
%% Visualization settings
fsize=12;   % Font size of figure
lw=2;       % linewidth of plot
%% Algorithm
for run = 1 : runs
    %% intialize weights
    U = zeros(1,order); % Input frame
    W = randn(1,order); % Initial Weigths
    %% Input/Output data
    N = 5000;               % Number of samples
    Bits = 2;               % Number of bits for modulation (2-bit for Binary modulation)
    data = randi([0 1],1,N);        % Random signal
    d = real(pskmod(data,Bits));    % BPSK Modulated signal (desired/output)
    r = filter(h,1,d);              % Signal after passing through channel
    x = awgn(r, SNRr);              % Noisy Signal after channel (given/input)
    for n = 1 : N
        U(1,2:end) = U(1,1:end-1);  % Sliding window
        U(1,1) = x(n);              % Present Input
     
        y = (W)*U';             % Calculating output of LMS
        e = d(n) - y;           % Instantaneous error 
        W = W +  eta * e * U ;  % Weight update rule of LMS
        J(run,n) = e * e';        % Instantaneous square error
    end
end
%% Calculation of performance parameters
MJ = mean(J,1);     % Mean square error
CS=freqz(h);        % Channel Spectrum
NF=(0:length(CS)-1)./(length(CS));          % Normalized Frequencies
IMR=-10*log10(real(CS).^2 + imag(CS).^2);   % Inverse channel magnitude response (desired)
IPR=-imag(CS)./real(CS);                    % Inverse channel phase response (desired)
ES=freqz(W);        % Equalizer Spectrum
EMR=10*log10(real(ES).^2 + imag(ES).^2);    % Equalizer magnitude response
EPR=imag(ES)./real(ES);                     % Equalizer phase response
%% Plots
figure % MSE
% plot(10*log10(MJ),'->k','linewidth',lw)
plot(10*log10(MJ),'-.k','linewidth',lw)
trendMJ = polyval(polyfit((0:N),[0 10*log10(MJ)],7),(1:N));
hold on
plot(trendMJ,'r','linewidth',lw)
hg=legend('MSE_{instantaneous}','MSE_{trend}','Location','Best');
grid minor
xlabel('Epochs iterations','FontSize',fsize);
ylabel('Mean squared error (dB)','FontSize',fsize);
title('Cost function','FontSize',2*fsize);
set(hg,'FontName','Times New Roman','FontSize',fsize)
set(gca,'FontName','Times New Roman','FontSize',fsize)
figure % Magnitude reponse
subplot(2,1,1)
plot(NF,IMR,'b','linewidth',lw)
hold on
plot(NF,EMR,'--r','linewidth',lw)
hg=legend('Inverse Channel','Equalizer','Location','Best');
grid minor
xlabel('Normalized Frequency','FontSize',fsize);
ylabel('Magnitude (dB)','FontSize',fsize);
title('Magnitude response','FontSize',2*fsize);
set(hg,'FontName','Times New Roman','FontSize',fsize)
set(gca,'FontName','Times New Roman','FontSize',fsize)
subplot(2,1,2)
plot(NF, IPR,'g','linewidth',lw)
hold on
plot(NF, EPR,'--b','linewidth',lw)
hg=legend('Inverse Channel','Equalizer','Location','Best');
grid minor
xlabel('Normalized Frequency','FontSize',fsize);
ylabel('Phase shift (rad)','FontSize',fsize);
title('Phase response','FontSize',2*fsize);
set(hg,'FontName','Times New Roman','FontSize',fsize)
set(gca,'FontName','Times New Roman','FontSize',fsize)