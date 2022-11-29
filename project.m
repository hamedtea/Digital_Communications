%Hamed Talebian - 150360360 - hamed.talebian@tuni.fi

%%% ---------------------------------------------------------------%%%
%%% the code and diagrams are optimized for different values of    %%%
%%% alpha(in Tx) and SNR. for not to observe plenty of unuseful    %%%
%%% diagrams, there are two methods:                               %%%
%%% 1- set the alpha/SNR parameter to a single value               %%%
%%% 2- run the code section by section                             %%%
%%% ---------------------------------------------------------------%%%
%%% The code is tested for 400000 bits.                            %%%
%%% to reduce time, required for compiling the code, change these: %%%
%%% N_symbols = 1000                                               %%%
%%% ---------------------------------------------------------------%%%

close all
clear all
clc

%general parameters
Rs = 40e6;           %symbol rate [sym/s]
T = 1/Rs;            %symbol time interval [s] 
r = 3;               %over-sampling factor
Fs = r/T;            % Sampling frequency 
Ts = 1/Fs;           % Sampling time interval
M = 16;              %16-QAM modulation size
%N_symbols = 1000;    %number of symbols
N_symbols = 100000;     % for relaxing the running time
alpha = 0.2:0.2:0.4;   %roll-off factor
Nd = 16;             %duration of RRC-filter in symbols
SNR = 0:1:20;      %signal to noise ratio
NFFT = 2^14;         %FFT size

%% Transmitter structure
%QAM alphabet generation
qam_axis = -sqrt(M)+1:2:sqrt(M)-1; %[-3 -1 1 3] 16-QAM axis
alphabet = bsxfun(@plus,qam_axis',1j*qam_axis);
alphabet = alphabet(:).';
alphabet_scaling_factor = 1/sqrt(mean(abs(alphabet).^2)); %mean power is 1
alphabet = alphabet*alphabet_scaling_factor;              %scaled alphabet generation

%QAM-Mapping 
N_bits = log2(M)*N_symbols; %number of bits
bits = randi(2,N_bits,1)-1; %column vector of N seudo-random bits

%% channel coding
k = 4; % number of source code
n = 7; %number of codewords
parity =[1 0 1; 1 1 1; 1 1 0; 0 1 1]; %define parity array - size(k*n-k)
%parity_array = parity_array;
%generator = [eye(k,k),parity_array];
%encoded_bits=zeros(1,n/k*length(bits).'); %matrix generation for encoded bits
id = eye(k);
generator=[parity id];
encoded_bits_final = []; %initialization of encoded bit sequence
encoded_bits = []; %initialization of final encoded bits
for i = 1:k:length(bits)-k+1
    bit_block = bits(i:i+k-1).';
    encoded_bits = mod(bit_block*generator,2); %encoding blocks of k bits into n coded bits
    encoded_bits_final = [encoded_bits_final; encoded_bits.']; %final encoded row vector %parity bits added at the begining of each block 
end

%% Gray coding
q = 2.^(log2(M)-1:-1:0); % [8 4 2 1]
B = reshape(encoded_bits_final,log2(M),[]); %blocks of 4 bits in columns 
symbol_indices = q*B;                       %bit blocks to decimal numbers - for gray coding 
[Gray_symbol_indices, mapgray] = bin2gray(symbol_indices, 'qam', M); 
symbols = alphabet(Gray_symbol_indices+1); %gray symbol generation


%% plotting the costellation  
figure (100)
plot(real(alphabet),imag(alphabet),'ro', 'MarkerFaceColor','r')
for i = 1:M
    text(real(alphabet(i))+ 0.02,imag(alphabet(i)), dec2base(mapgray(i),2,4)) %0.02 for clarity of plot
end
axis equal
grid on
xlabel(' I (Re)') 
ylabel('Q (Im)') 
title('Gray coded 16-QAM constellation')
%% Upsampling and RRC pulse shaping for differet values of alpha
symbols_upsampled = zeros(1, length(r*symbols)); % Zero vector for up-sampled symbol sequence 
symbols_upsampled(1:r:r*length(symbols)) = symbols;  % plugging samples into unsampled vector [a0 0 0 a1 0 0  ...]

% FIR RRC filter generation
for i= 1:1:length(alpha)
    pulse(i,:) = rcosdesign(alpha(i), Nd, r, 'sqrt'); %nyquist pulse shape %alpha = roll-off factor Nd = number of symbols per pulse r = oversampling factor 'sqrt = normal RRC
    st_temp = filter(pulse(i,:),1,symbols_upsampled); % Signal transmitted after filtering : pulse shaping 
    filter_delay = (length(pulse(i,:))-1)/2;
    st(i,:) = st_temp(1+filter_delay:end);     % Filter delay correction
end
%% plotting Tx eye diagram for different alpha values
for i= 1:1:length(alpha)
    figure (1+i)
    hold on
    for ii = 1:2*r:(length(st(i,1:1:2000)))  %first 5000 samples is selected for plotting %imaginary part is the same (not to plot)
        plot(real(st(i,ii:ii+2*r)),'b') 
        %plot real part of Tx signal for each upsamled symbol duration 
    end 
    hold off 
    grid on
    xlabel('r [T/Ts]');ylabel('Amplitude')
    title(['Tx eye diagram (real part), \alpha = ', num2str(alpha(1,i))])
end
%% plotting to probe the role of exess bandwidth - for fun
f = -Fs/2:1/(NFFT*Ts):Fs/2-1/(NFFT*Ts); %frequency vector for plotting purposes
for i= 1:1:length(alpha)
    figure(999+i)
    subplot(221)
    plot((-Nd*r*Ts:Ts:Nd*r*Ts),conv(pulse(i,:),pulse(i,:)))
    xlabel('Time [Ts]');ylabel('Relative Amplitude')
    title(['Tx/Rx RRC convoluton, \alpha = ', num2str(alpha(i)), ''])
    legend('Pulse shape','Ideal symbol-sampling locations')
    subplot(222)
    impz(pulse(i,:))
    title(['Tx/Rx RRC impulse response, \alpha = ', num2str(alpha(i)), ''])
    xlabel('samples');ylabel('relative Amplitude')
    subplot(223)
    plot((-Nd*r/2*Ts:Ts:Nd*r/2*Ts),pulse(i,:));hold on;grid on
    stem(-Nd*r/2*Ts:T:Nd*r/2*Ts,pulse(i,1:r:end),'ro') %sampling times
    xlabel('Time [Ts]');ylabel('Relative Amplitude')
    title(['RRC nyquist pulse shape, \alpha = ', num2str(alpha(i)), ''])
    legend('Pulse shape','Ideal symbol-sampling locations')
    hold off
    subplot(224)
    plot(f, fftshift(abs(fft(pulse(i,:),NFFT))));
    xlabel('Frequency [MHz]');ylabel('Relative Amplitude ')
    ylabel('amplitude [db]')
    title(['RRC nyquist pulse in freq. domain, \alpha = ', num2str(alpha(i)), '']);hold on;grid on
end
%% plotting spectral contents of transmitted signal 
%and plotting RRC pulse in time/frequency domain
f = -Fs/2:1/(NFFT*Ts):Fs/2-1/(NFFT*Ts); %frequency vector for plotting purposes
for i= 1:1:length(alpha)
    figure (100+i)
    subplot(2,2,1) 
    plot(f/1e6, fftshift(abs(fft(st(i,:), NFFT)))); %spectrum of recieved signal
    grid on
    xlabel('Frequency [MHz]') 
    ylabel('Relative Amplitude') 
    title(['TX signal, after Tx filtering, \alpha = ', num2str(alpha(i))]) 
    subplot(2,2,2); 
    plot(f/1e6, fftshift(20*log10(abs(fft(st(i,:), NFFT))))); %spectrum of recieved signal [db format]
    grid on
    xlabel('Frequency [MHz]') 
    ylabel('Amplitude [db] ') 
    title(['TX signal, after Tx filtering \alpha = ', num2str(alpha(i))])
    subplot(2,2,3) 
    %impz(pulse)
    plot(-Nd*r/2*Ts:Ts:Nd*r/2*Ts,pulse(i,:),'b') %RRC pulse 
    hold on
    grid on
    stem(-Nd*r/2*Ts:T:Nd*r/2*Ts,pulse(i,1:r:end),'ro') %sampling times
    xlabel('time [s]') 
    ylabel('Amplitude')
    title(['Transmit/receive RRC filter (pulse shape), \alpha = ', num2str(alpha(i))])
    legend('Pulse shape','Ideal symbol-sampling locations')
    subplot(2,2,4)
    plot(f/1e6, fftshift(10*log10(abs(fft(pulse(i,:), NFFT))))); %amplitude response of RRC / RRC in freq. domain
    xlabel('Frequency [MHz]') 
    ylabel('Amplitude [db]') 
    grid on
    title(['Amplitude response of RRC filter, \alpha = ', num2str(alpha(i))])
    hold off
end

%% Channel modeling
%multipath channel model filtering
h = [0.19+.56j .45-1.28j -.14-.53j -.19+.23j .33+.51j]; % Complex channel model (coefficients)
strong_tap = abs(h); %finding the strongest tap ---> second one
L1=1; L2=3; % Channel maximum tap is the second one
[H_ch,f_ch]=freqz(h,1,-Rs/2:Rs/400:Rs/2,Rs); %channel freq/(angle) phase response

%% %Plotting amplitude and phase response of the channel: 
figure(2); 
subplot(311); 
stem(-L1:L2,abs(h),'r'); % Absolute values of channel impulse response / FIR of equivalent channel model
xlabel('Time delay');ylabel('Absolute Coefficients');
title('Absolute values of impulse responses'); 
subplot (312); 
plot(f_ch/1e6,20*log10(abs(H_ch)),'b'); % FIR freq response 
grid on; 
xlabel('Frequency [MHz]');ylabel('Amplitude response [dB]'); 
title('Frequency Responses');
subplot (313);
plot(f_ch/1e6,angle(H_ch)*180/pi,'r'); % phase response of equivalent channel model
grid on; 
xlabel('Frequency [MHz]');
ylabel('Angle [degree]'); 
title('Phase Responses'); 

%% AWGN model

mt = filter(h,1,st(1,:)); % signal after ISI FIR filtering 
mt = mt(L1:end); % discard the delay
%AWGN noise addition
    noise = (1/sqrt(2))*(randn(size(mt)) + 1j*randn(size(mt))); %random complex noise generator
    P_s = var(mt);              % Signal power 
    P_n = var(noise);           % Noise power   
    noise_scaling_factor = sqrt(P_s/P_n./10.^(SNR./10)*(r/(1+alpha(1,1))));

%define recieved signal for each SNR and alpha value
mt_noisy = zeros(length(SNR), length(mt)); %3D matrix: (:alpha,:SNR,:noisy_samples)

for ii=1:1:length(SNR)
    % Defining noise scaling factor based on the desired SNR: 
    mt_noisy(ii,:) = mt + noise_scaling_factor(ii)*noise; %recieved signal after adding each AWGN noise
end

%% plotting recieved signal spectrum after channel effects (adding IS and noise) for SNR values
f = -Fs/2:1/(NFFT*Ts):Fs/2-1/(NFFT*Ts); %frequency vector for plotting purposes
for ii=1:10:length(SNR)
    figure(200+ii)
    plot(f/1e6, fftshift(20*log10(abs(fft(mt_noisy(ii,:), NFFT))))); %spectral content of recieved signal
    hold on
    plot(f_ch/1e6, 20*log10(H_ch), 'r') %corresponding fading channel frequency response
    xlabel('Frequency [MHz]') 
    ylabel('Amplitude [db] ') 
    legend('recived signal spectrum', 'Channel Amplitude response')
    title(['signal spectrum plus channel effects, SNR =', num2str(SNR(ii)),''])
    hold off
end
%% Reciever structure - signal with noise and multipaths
%Rx filter
alpha_index = find(alpha==0.2); %find the index of alpha = 0.2 for setting Rx filter
qtt = zeros(length(SNR), length(mt_noisy)-filter_delay);    %2D matrix of reieved samples for each value of alpha and SNR
%filtering and down sampling for each alpha & SNR values
for ii=1:1:length(SNR)
    ft = pulse(alpha_index,:);             % signal at the begining of reciever : suitable for probing the channel effects
    qt = filter(ft,1,mt_noisy(ii,:));      % Receiver filtering 
    qt = qt(1+filter_delay:end);           % Filter delay correction
    qtt(ii,:) = qt;                        % adding each alpha/SNR filtered signal into a 3D matrix
    qk(ii,:)=qtt(ii,1:r:end);              % down sampling for each alpha/SNR values a 3D matrix for equalization
end

%% Plotting and filtering the transmitted signal without channel effects (multipaths/noise) to visualize ideal eye diagram
qk_Not_effected = zeros(length(alpha),length(symbols)-Nd); %3D matrix for each value of alpha and SNR
qtt_Not_effected = zeros(length(alpha),length(st)-filter_delay);
for i= 1:1:length(alpha)
    ft = pulse(i,:); % signal at the begining of reciever : suitable for probing the channel effects
    qt = filter(ft,1,st(i,:));      % Receiver filtering 
    qt = qt(1+filter_delay:end); % Filter delay correction
    qtt_Not_effected(i,:) = qt;
    %qk_Not_effected(i,:)=qtt_Not_effected(i,:); %down sampling for SNR values
end

for i= 1:1:length(alpha)
    figure (300+i)
    subplot(211)
    hold on 
    for iii = 1:2*r:(length(qtt_Not_effected(i,1:1:1000)))
        plot(real(qtt_Not_effected(i,iii:iii+2*r)),'b'); 
    end 
    hold off
    xlabel('time')
    ylabel('Amplitude')
    grid on
    title(['Rx eye diagram (without channel effects) - real part, \alpha = ', num2str(alpha(i)), ''])
    subplot(212)
    hold on
    for iii = 1:2*r:(length(qtt_Not_effected(i,1:1:1000)))
        plot(imag(qtt_Not_effected(i,iii:iii+2*r)),'b');  
    end 
    hold off
    grid on
    xlabel('time')
    ylabel('Amplitude')
    grid on
    title(['Rx eye diagram (without channel effects) - imaginary part, \alpha = ', num2str(alpha(i)), ''])
end
%% Plotting Rx signal eye diagram for SNR values
for ii=1:10:length(SNR)
    figure (400+ii)
    subplot(211)
    hold on 
    for iii = 1:2*r:(length(qtt(ii,1:1:1000))) 
            plot(real(qtt(ii,iii:iii+2*r)),'b'); 
    end 
    hold off
    xlabel('time')
    ylabel('Amplitude')
    grid on
    title(['Rx eye diagram - real part, SNR = ', num2str(SNR(:,ii)), ''])
    subplot(212)
    hold on
    for iii = 1:2*r:(length(qtt(ii,1:1:1000))) 
            plot(imag(qtt(ii,iii:iii+2*r)),'b');  
    end 
    hold off
    grid on
    xlabel('time')
    ylabel('Amplitude')
    grid on
    title(['Rx eye diagram - imaginary part, SNR = ', ...
            num2str(SNR(:,ii)), ''])
end  
%% plotting recieved constellation for SNR values
for ii=1:10:length(SNR)
   figure(500+ii)
   plot(real(qk(ii,:)), imag(qk(ii,:)), 'b*');
   hold on; 
   grid on;
   plot(symbols,'ro', 'MarkerFaceColor','r'); 
   legend('Received samples', 'Original symbols')
   xlabel('Real Part')
   ylabel('Imaginary Part')
   title(['Received symbol-rate samples (RX constellation), ' ...
            'SNR = ', num2str(SNR(ii))])
   axis equal
   hold off
end

%% LMS Equalizer 
mu = 0.01; % step-size of the algorithm 
order = 20;
iteration = 10000;
%c_LMS = zeros(order+1,1); % equalizer coefficients, initialization
C_LMSS = zeros(length(SNR),order+1); % equalizer coefficients, initializations 
EKK = [];
for ii=1:1:length(SNR)
    disp 'SNR = ', disp(SNR(ii))
    c_LMS = zeros(order+1,1);
    for iii = iteration+1:length(qk)-iteration
        rk = flipud(qk(ii,iii-order/2:iii+order/2).'); % Received signal vector :transpose/flipud
        e(iii) = symbols(iii) - c_LMS.'*rk; % Error signal, we assume a known symbol sequence 
        c_LMS = c_LMS + mu*e(iii)*conj(rk); % LMS update !
    end
    C_LMSS(ii,:) = c_LMS;
    EKK(ii,:) = e;
end
%% plotting convergence behaviour of LMS for SNR values
for  ii=1:5:length(SNR)
     figure(600+ii)
     plot(abs(EKK(ii,1000:3000)))
     ylabel('LMS error')
     xlabel('Iteration index')
     title(['Convergence behavior of the LMS-algorithm, SNR =', num2str(SNR(ii)), '']);
     hold off
end

%% plotting effective impulse response of equalized system for different SNR and alpha values
for ii=1:10:length(SNR)
    figure(700+ii) 
    hold on 
    stem(abs(conv(h,C_LMSS(ii,:)))); %plot impulse response LMS equilezer
    title(['Effective impulse response (abs) of the equalized system, SNR = ', num2str(SNR(:,ii)), '']) 
    hold off
end

%% plotting and comparing channel frequency and phase response of the LMS equalizer
for ii=1:10:length(SNR)
    [H_c_LMS,f_c_LMS]=freqz(C_LMSS(ii,:),1,-Rs/2:Rs/400:Rs/2,Rs);
    [H_tot,f_tot]=freqz(conv(C_LMSS(ii,:),h),1,-Rs/2:Rs/400:Rs/2,Rs); 
    figure(800+ii)
    plot(f_ch/1e6,20*log10(abs(H_ch)),'b')
    hold on
    plot(f_c_LMS/1e6,20*log10(abs(H_c_LMS)),'r');
    plot(f_tot/1e6,20*log10(abs(H_tot)),'g');
    legend('Channel','LMS Equalizer','Total Response (LMS)')
    xlabel('Freuency [MHz]')
    ylabel('Amplitude [dB]')
    title(['LMS Equalizer performance (frequency responce) SNR = ', ...
             num2str(SNR(ii)), ''])
    hold off
 end


%% Plotting for comparing phase response the LMS equalizer
for ii=1:10:length(SNR)
    figure (900+ii)
    plot(f_ch/1e6,angle(H_ch)*180/pi,'b')
    hold on
    plot(f_ch/1e6,angle(H_c_LMS)*180/pi,'r')
    plot(f_ch/1e6,angle(H_tot)*180/pi,'g')
    grid on 
    xlabel('Frequency [MHz]');
    ylabel('Angle [degree]'); 
    legend('Channel','LMS Equalizer','Total Response (LMS)'); 
    title(['LMS Equalizer performance (Phase responce) SNR = ', ...
    num2str(SNR(ii)), '']);
end
%% filtering with equalizer coefficients
qk_equalized_final = zeros(length(SNR), length(qk)-order); 
                        %-length(C_LMSS(:,length(SNR))),length(SNR), length(alpha)); % columns = length qk - delay qk (lenght CMS) & rows = SNR values
for ii=1:1:length(SNR)
    qk_equalized = filter(C_LMSS(ii,:),1,qk(ii,:));
    qk_equalized = qk_equalized(1+order:end); % Filter delay correction order--> order of equalizerfor delay correction
    qk_equalized_final(ii,:)= qk_equalized;
end
%% plotting equalized constellation - for 5 db SNR difference 
for ii=1:5:length(SNR)
    figure(1000+ii) 
    plot(real(qk(ii,1:10000)), imag(qk(ii,1:10000)),'*b') %select the proper qk-array - squeeze tranpose to back to column vector 
    hold on 
    plot(real(qk_equalized_final(ii,1:10000)), imag(qk_equalized_final(ii,1:10000)), 'go') 
    plot(real(symbols), imag(symbols), 'ro', 'MarkerFaceColor','r' );
    legend('Received constellation','LMS Equalized constellation', 'original symbols')
    xlabel('Real part')
    ylabel('Imaginary part')
    title(['Equalized constellation, SNR = ', ...
    num2str(SNR(ii)), ''])
    hold off
end
%% Symbol-by-symbol minimum distance detection and gray decoding for each SNR value
estimated_bits = [];
for ii=1:1:length(SNR)
    alphabet_error_matrix = abs(bsxfun(@minus,alphabet.',qk_equalized_final(ii,:))); 
        %rows = alphabet & columns = recieved equalized symbol indices
    [~,estimated_gray_symbol_ind] = min(alphabet_error_matrix); %search for minimum distance
    %Gray decoding
    estimated_symbol_indices = gray2bin(estimated_gray_symbol_ind-1,'qam',M);
    estimated_bit_blocks = rem(floor((estimated_symbol_indices(:))*2.^(1-log2(M):0)),2).'; %rows = SNR columns bit by bit minimum distance detection
    estimated_bits(ii,:) = estimated_bit_blocks(:);
end

%% channel decoding
p = n-k; %size of identity matrix
parity_check = [eye(p) parity.']; %parity check matrix  
parity_check= parity_check.'; %transpose parity check 
decoded_bits= []; %column array generation for decoded bits
s=[]; %syndrome initialization
decoded_bits_final = [];
[~, numofcolumns] = size(estimated_bits); %retrive the length for loops
for ii=1:1:length(SNR)
    for j=1:n:numofcolumns-n+1
        codewords = estimated_bits(ii,j:j+n-1);
        s = mod(codewords(1,:)*parity_check,2);
        if s == 0
            valid_code=codewords((n-k+1):n);
        else
            pos=0; %initialization of error position finder for hard decoding
            for mm=1:length(parity_check)
            if(s(1,:)==parity_check(mm,:)) %finding the poisition of error
                pos=mm;
            end
            end
        zero=zeros(1,n);
        zero(pos)=~zero(pos);
        errorpat=zero; %error pattern variable
        error_correction=xor(codewords,errorpat); %modulu correction
        valid_code = error_correction((n-k+1):n);
        end
        decoded_bits = [decoded_bits;valid_code.'];
    end
    decoded_bits_final(ii,:) = decoded_bits; %final decoded matrix - rows: SNR values columns:decoded/corrected bits
    decoded_bits = []; %empthy row vector for next SNR iteration 
end

%% BER calculations

bit_errors = decoded_bits_final ~= bits(1:length(decoded_bits_final)).'; %finding errors 
BER = mean(bit_errors, 2); %mean values of bit errors 

figure (111)
semilogy(SNR, BER, 'LineWidth', 3);
hold on; 
title('Bit error rate') 
xlabel('SNR [dB]') 
ylabel('BER [log scale]') 
legend('Simulated BER','Theoretical symbol error probability');

