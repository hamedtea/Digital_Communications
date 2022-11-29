%Hamed Talebian - 150360360 - hamed.talebian@tuni.fi

%%% ---------------------------------------------------------------%%%
%%% the code and diagrams are optimized for different values of    %%%
%%% alpha and SNR. for not to observe plenty of unuseful diagrams, %%%
%%% there are two methods:                                         %%%
%%% 1- set the alpha/SNR parameter to a single value               %%%
%%% 2- run the code section by section                             %%%
%%% ---------------------------------------------------------------%%%
%%% The code is tested for 400000 bits.                            %%%
%%% to reduce time, required for compiling the code, change these: %%%
%%% N_symbols = 1000 
%%% iteration = 30 %to equalize less number of symbols
%%% order = 12
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
N_symbols = 1000;    %number of symbols
%N_symbols = 100000;     % for relaxing the running time
alpha = 0.2;   %roll-off factor
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
bit_blocks = reshape(bits,[], k); %blocks of 4 bits to add parity check / part of Gray coding is also done in here
encoded_bits = []; %initialization of encoded bit sequence
encoded_bits_final = [];
for i=1:length(bit_blocks)
    encoded_bits = [mod(bit_blocks(i,:)*generator,2)]; %encoding blocks of k bits into n coded bits
    t = reshape(encoded_bits, [],1);
    encoded_bits_final = [encoded_bits_final; t];
end

%% Gray coding
q = 2.^(log2(M)-1:-1:0); % [8 4 2 1]
B = reshape(encoded_bits_final,log2(M),[]); %it is done in channel coding part 
symbol_indices = q*B;
[Gray_symbol_indices, mapgray] = bin2gray(symbol_indices, 'qam', M); 
[Gray_symbol_indices, mapgray] = bin2gray(symbol_indices, 'qam', M); 
symbols = alphabet(Gray_symbol_indices+1); %gray symbol generation


%% plotting the costellation  
figure (100)
plot(real(alphabet),imag(alphabet),'ro')
for i = 1:M
    text(real(alphabet(i)),imag(alphabet(i)), dec2base(mapgray(i),2,4))
end
axis equal 
xlabel(' I (Re)') 
ylabel('Q (Im)') 
title('Gray coded 16-QAM constellation')
%% Upsampling and RRC pulse shaping for differet values of alpha

symbols_upsampled = zeros(1, length(r*symbols)); % Zero vector for up-sampled symbol sequence 
symbols_upsampled(1:r:r*length(symbols)) = symbols;  % plugging samples into unsampled vector [a0 0 0 a1 0 0  ...]

% FIR RRC filter generation
gt = rcosdesign(alpha, Nd, r, 'sqrt'); %nyquist pulse shape %alpha = roll-off factor Nd = number of symbols per pulse r = oversampling factor 'sqrt = normal RRC


%Symbol insertation
st = filter(gt,1,symbols_upsampled);
filter_delay = (length(gt)-1)/2;
ft = ft(1+filter_delay:end);     % Filter delay correction

%% plotting Tx eye diagram for different alpha values
for ii = 1:2*r:(length(st(i,:))-2*r)  
    figure(2)
    plot(real(st(i,ii:ii+2*r)),'b') 
    %plot real part of Tx signal for each upsamled symbol duration 
    hold on
end 
hold off 
grid on
xlabel('r [T/Ts]');ylabel('Amplitude')
title(['Tx eye diagram (real part), \alpha = ', num2str(alpha(1,i))])



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
for i= 1:1:length(alpha) %filter the transmitted signals for alpha values %columns = pulses rows = each alpha value
    mt = filter(h,1,st(i,:)); % signal after ISI FIR filtering 
    mt = mt(L1:end); % discard the delay
%AWGN noise addition
    noise = (1/sqrt(2))*(randn(size(mt)) + 1j*randn(size(mt))); %random complex noise generator
    P_s = var(mt);              % Signal power 
    P_n = var(noise);           % Noise power   
    noise_scaling_factor(i,:) = sqrt(P_s/P_n./10.^(SNR./10)*(r/(1+alpha(i))));
end

%define recieved signal for each SNR and alpha value
mt_noisy = zeros(length(alpha), length(SNR), length(mt)); %3D matrix: (:alpha,:SNR,:noisy_samples)
for i= 1:1:length(alpha)
    for ii=1:1:length(SNR)
    % Defining noise scaling factor based on the desired SNR: 
    mt_noisy(i,ii,:) = mt + noise_scaling_factor(i,ii)*noise; %recieved signal after adding each AWGN noise
    end
end

%% plotting recieved signal spectrum after channel effects (adding IS and noise) for SNR values
f = -Fs/2:1/(NFFT*Ts):Fs/2-1/(NFFT*Ts); %frequency vector for plotting purposes
for i= 1:1:length(alpha)
    for ii=1:5:length(SNR)
        sq_mt_noisy = squeeze(mt_noisy(i,ii,:)).'; %selecting the 3 dimention, corresponding to each alpha value
        figure(200+j+ii)
        plot(f/1e6, fftshift(20*log10(abs(fft(sq_mt_noisy, NFFT))))); %spectral content of recieved signal
        hold on
        plot(f_ch/1e6, 20*log10(H_ch), 'r') %corresponding fading channel frequency response
        xlabel('Frequency [MHz]') 
        ylabel('Amplitude [db] ') 
        legend('recived signal spectrum', 'Channel Amplitude response')
        title(['signal spectrum plus channel effects, SNR =', ...
            num2str(SNR(ii)), '\alpha = ' num2str(alpha(i))])
        hold off
    end
    j=j+length(SNR); %for plotting purpos-a new figure in each iteration
end
%% Reciever structure - signal with noise and multipaths
%Rx filter
qk = zeros(length(alpha), length(SNR), length(symbols)-Nd);                %3D matrix for each value of alpha and SNR
qtt = zeros(length(alpha), length(SNR), length(mt_noisy)-filter_delay);    %3D matrix of reieved samples for each value of alpha and SNR
%filtering and down sampling for each alpha & SNR values
for i= 1:1:length(alpha)
    for ii=1:1:length(SNR)
        ft = pulse(i,:);                         % signal at the begining of reciever : suitable for probing the channel effects
        %mt_s = squeeze(mt_noisy(i,ii,:)).'; %squeezeto a column vector for each SNR/ALpha value
        qt = filter(ft,1,mt_noisy(i,ii,:));      % Receiver filtering 
        qt = qt(1+filter_delay:end);             % Filter delay correction
        qtt(i,ii,:) = qt;                        % adding each alpha/SNR filtered signal into a 3D matrix
        qk(i,ii,:)=qtt(i,ii,1:r:end);            %down sampling for each alpha/SNR values a 3D matrix for equalization
    end
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
    for iii = 1:2*r:(length(qtt_Not_effected(i,:))-2*r) 
        plot(real(qtt_Not_effected(i,iii:iii+2*r)),'b'); 
    end 
    hold off
    xlabel('time')
    ylabel('Amplitude')
    grid on
    title(['Rx eye diagram (without channel effects) - real part, \alpha = ', num2str(alpha(i)), ''])
    subplot(212)
    hold on
    for iii = 1:2*r:(length(qt)-2*r) 
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
for i= 1:1:length(alpha)
    for ii=1:5:length(SNR)
        figure (400+j+ii)
        subplot(211)
        hold on 
        for iii = 1:2*r:(length(qtt(i,ii,:))-2*r) 
            plot(squeeze(real(qtt(i,ii,iii:iii+2*r))).','b'); 
        end 
        hold off
        xlabel('time')
        ylabel('Amplitude')
        grid on
        title(['Rx eye diagram - real part, SNR = ', num2str(SNR(:,ii)), ...
            ', \alpha = ', num2str(alpha(i)), ''])
        subplot(212)
        hold on
        for iii = 1:2*r:(length(qt)-2*r) 
            plot(squeeze(imag(qtt(i,ii,iii:iii+2*r))).','b');  
        end 
        hold off
        grid on
        xlabel('time')
        ylabel('Amplitude')
        grid on
        title(['Rx eye diagram - imaginary part, SNR = ', ...
            num2str(SNR(:,ii)), ', \alpha = ', num2str(alpha(i)), ''])
    end
    j=j+length(SNR); %dummy variable for plotting purpose
end

%% plotting recieved constellation for SNR/alpha values
for i= 1:1:length(alpha)
    for ii=1:5:length(SNR)
        figure(500+j+ii)
        plot(real(squeeze(qk(i,ii,:)).'), imag(squeeze(qk(i,ii,:)).'), 'b*');
        %plot(qk(i,ii,:),'b*'); 
        hold on; 
        grid on;
        plot(symbols,'ro', 'MarkerFaceColor','r'); 
        legend('Received samples', 'Original symbols')
        xlabel('Real Part')
        ylabel('Imaginary Part')
        title(['Received symbol-rate samples (RX constellation), ' ...
            'SNR = ', num2str(SNR(ii)), ' \alpha = ', num2str(alpha(i)) ])
        axis equal
        hold off
    end
    j=j+length(SNR); %dummy variable for plotting purpose
end

%% LMS Equalizer 
mu = 0.01; % step-size of the algorithm 
order = 30;
iteration = 60;
%c_LMS = zeros(order+1,1); % equalizer coefficients, initialization
C_LMSS = zeros(length(alpha),length(SNR),order+1); % equalizer coefficients, initializations 
EKK = zeros(length(alpha), length(SNR),length(qk)-order);
for i= 1:1:length(alpha)
    for ii=1:1:length(SNR)
        disp 'SNR = ', disp(SNR(ii))
        c_LMS = zeros(order+1,1);
        qk_s = squeeze(qk(i,ii,:)).';
        for iii = iteration+1:length(qk_s)-iteration
            rk = flipud(qk_s(iii-order/2:iii+order/2).'); % Received signal vector :transpose/flipud
            EK(iii) = symbols(iii) - c_LMS.'*rk; % Error signal, we assume a known symbol sequence 
            c_LMS = c_LMS + mu*EK(iii)*conj(rk); % LMS update !
        end
    C_LMSS(i,ii,:) = c_LMS;
    EKK(i,ii,:) = EK;
    end
end

%% plotting convergence behaviour of LMS for SNR values
for  ii=1:5:length(SNR)
    figure(600+i+j)
    for i= 1:1:length(alpha)
        plot(abs(squeeze(EKK(i,ii,:)).'))
        ylabel('LMS error')
        xlabel('Iteration index')
        hold on
    end
    Legend=cell(length(alpha),1);
    for iter=1:length(alpha)
        Legend{iter}=strcat('\alpha =', num2str(alpha(iter)));
    end
    legend(Legend)
    title(['Convergence behavior of the LMS-algorithm, SNR =', num2str(SNR(ii)), '']);
    hold off
    j=j+length(SNR); %dummy variable for plotting purpose
end

%% plotting effective impulse response of equalized system for different SNR and alpha values
for ii=1:5:length(SNR)
    figure(700+ii+j) 
    for i= 1:1:length(alpha)
        hold on 
        squeeze_C_LMSS = squeeze(C_LMSS(i,ii,:)).'; %reduce 3D matrix dimention for plotting 
        stem(abs(conv(h,squeeze_C_LMSS))); %plot impulse response LMS equilezer
        title(['Effective impulse response (abs) of the equalized system, SNR = ', num2str(SNR(:,ii)), '']) 
        hold off
    end
    Legend=cell(length(alpha),1);
    for iter=1:length(alpha)
        Legend{iter}=strcat('\alpha =', num2str(alpha(iter)));
    end
    legend(Legend)
    j=j+length(SNR); %dummy variable for plotting purpose
end

%% plotting and comparing channel frequency and phase response of the LMS equalizer
for i= 1:1:length(alpha)     
    for ii=1:5:length(SNR)
        squeeze_C_LMSS = squeeze(C_LMSS(i,ii,:)).'; %%reduce 3D matrix dimention for plotting
        [H_c_LMS,f_c_LMS]=freqz(squeeze_C_LMSS,1,-Rs/2:Rs/400:Rs/2,Rs);
        [H_tot,f_tot]=freqz(conv(squeeze_C_LMSS,h),1,-Rs/2:Rs/400:Rs/2,Rs); 
        figure(800+j+ii)
        plot(f_ch/1e6,20*log10(abs(H_ch)),'b')
        hold on
        plot(f_c_LMS/1e6,20*log10(abs(H_c_LMS)),'r');
        plot(f_tot/1e6,20*log10(abs(H_tot)),'g');
        legend('Channel','LMS Equalizer','Total Response (LMS)')
        xlabel('Freuency [MHz]')
        ylabel('Amplitude [dB]')
        title(['LMS Equalizer performance (frequency responce) SNR = ', ...
             num2str(SNR(ii)), ' \alpha = ', num2str(alpha(i)) ])
        hold off
    end
    j=j+length(SNR); %dummy variable for plotting purpose
end

%% Plotting for comparing phase response the LMS equalizer
for i= 1:1:length(alpha) 
    for ii=1:5:length(SNR)
        figure (900+j+ii)
        plot(f_ch/1e6,angle(H_ch)*180/pi,'b')
        hold on
        plot(f_ch/1e6,angle(H_c_LMS)*180/pi,'r')
        plot(f_ch/1e6,angle(H_tot)*180/pi,'g')
        grid on 
        xlabel('Frequency [MHz]');
        ylabel('Angle [degree]'); 
        legend('Channel','LMS Equalizer','Total Response (LMS)'); 
        title(['LMS Equalizer performance (Phase responce) SNR = ', ...
            num2str(SNR(ii)), ' \alpha = ', num2str(alpha(i)) ]);
    end
    j=j+length(SNR); %dummy variable for plotting purpose
end
%% filtering with equalizer coefficients
qk_equalized_final = zeros(length(alpha), length(SNR), length(qk)-order); 
                        %-length(C_LMSS(:,length(SNR))),length(SNR), length(alpha)); % columns = length qk - delay qk (lenght CMS) & rows = SNR values
for i= 1:1:length(alpha)
    for ii=1:1:length(SNR)
        sq_C_LMSS=squeeze(C_LMSS(i,ii,:)).'; %reduce dimenstion for filtering
        sq_qk = squeeze(qk(i,ii,:)).'; %reduce dimenstion for filtering
        qk_equalized = filter(sq_C_LMSS,1,sq_qk);
        qk_equalized = qk_equalized(1+order:end); % Filter delay correction order--> order of equalizerfor delay correction
        qk_equalized_final(i,ii,:)= qk_equalized;
    end
end
%% plotting equalized constellation - for 5 db SNR difference 
for i= 1:1:length(alpha)
    for ii=1:5:length(SNR)
        figure(1000+j+ii) 
        sq_qk_equalized_final=squeeze(qk_equalized_final(i,ii,:)).'; %reduce dimenstion for filtering
        sq_qk = squeeze(qk(i,ii,:))'; %reduce dimenstion for filtering
        plot(real(sq_qk), imag(sq_qk),'*b') %select the proper qk-array - squeeze tranpose to back to column vector 
        hold on 
        plot(real(sq_qk_equalized_final), imag(sq_qk_equalized_final), 'go') 
        plot(real(symbols), imag(symbols), 'ro', 'MarkerFaceColor','r' );
        legend('Received constellation','LMS Equalized constellation', 'original symbols')
        xlabel('Real part')
        ylabel('Imaginary part')
        title(['Equalized constellation, SNR = ', ...
            num2str(SNR(ii)), ' \alpha = ', num2str(alpha(i))])
        hold off
    end
    j=j+length(SNR); %dummy variable for plotting purpose
end
%% Symbol-by-symbol minimum distance detection and gray decoding for each SNR value
estimated_bits = [];
for i= 1:1:length(alpha)
    for ii=1:1:length(SNR)
        alphabet_error_matrix = abs(bsxfun(@minus,alphabet.',qk_equalized_final(i,ii,:))); 
        %rows = alphabet & columns = recieved equalized symbol indices
        [~,estimated_gray_symbol_ind] = min(alphabet_error_matrix); %search for minimum distance
    %Gray decoding
        estimated_symbol_indices = gray2bin(estimated_gray_symbol_ind-1,'qam',M);
        estimated_bit_blocks = rem(floor((estimated_symbol_indices(:))*2.^(1-log2(M):0)),2).'; %rows = SNR columns bit by bit minimum distance detection
        estimated_bits(i,ii,:) = estimated_bit_blocks(:);
        %estimated_coded_bits(:,ii,i) = estimated_bits; %rows = encoded_bits columns = SNR values
    end
end

%% channel decoding
p = n-k; %size of identity matrix
parity_check = [eye(p) parity.']; %parity check matrix  
parity_check= parity_check.'; %transpose parity check 
decoded_bits= []; %column array generation for decoded bits
s=[]; %syndrome initialization
decoded_bits_final = [];
valid_code_temp = [];
%for i= 1:1:length(alpha)
%dcbits = reshape(estimated_bits,[],n);
for i=1:1:length(alpha)
for ii=1:1:length(SNR)
    decoded_temp = estimated_bits(i,ii,:); %select a column array from 3D matrix for parity check
    for jj=1:n:length(decoded_temp)-n
        r = decoded_temp(i,ii,jj:jj+n-1);
        s = mod(r(i,ii,:).'*parity_check,2);
        if s == 0
            %validcode_f=r((n-k+1):n);
            validcode_f=r(1:k);
        else
            pos=0; %initialization of error position finder for hard decoding
            for mm=1:length(parity_check)
            if(s(1,:)==parity_check(mm,:)) %finding the poisition of error
                pos=mm;
            end
            end
        zero=zeros(1,n);
        zero(pos)=~zero(pos);
        errorpat=zero;
        validcode=xor(r,errorpat);
        %valid_code_f = validcode((n-k+1):n);
        valid_code_f = validcode(1:k);
        end
        tt = reshape(validcode_f,[],1);
        decoded_bits = [decoded_bits;tt];
    end
    %decoded_bits_final(ii,:) = reshape(decoded_bits,[],ii);
    decoded_bits_final(i,ii,:) = decoded_bits(:); %final column array
    decoded_bits = []; 
end
end
%decoded_bits_final = reshape(decoded_bits, length(SNR), []);
%end
%decoded_bits_final = reshape(decoded_bits, length(SNR), []); %reshape decoded bits (column vector) based on SNR values

%% BER calculations
d = min(abs(alphabet(1)-alphabet(2:end)));
sigma = sqrt(0.5 * P_n * noise_scaling_factor.^2); 
P_sym_error = (4*qfunc(d./(2*sigma)) - 4*qfunc(d./(2*sigma)).^2) ... 
    * ((M-4-4*(sqrt(M)-2))/M) + ... 
    (2*qfunc(d./(2*sigma))- qfunc(d./(2*sigma)).^2) * (4/M) + ... 
    (3*qfunc(d./(2*sigma)) - 2*qfunc(d./(2*sigma)).^2) * (4*(sqrt(M)-2)/M); 

BER = zeros(1,length(SNR));
for ii=1:1:length(SNR)
    %bit_errors = estimated_coded_bits(:,ii) ~= bits(1:length(estimated_coded_bits)); 
     bit_errors = decoded_bits_final(1,ii,:) ~= bits(length(decoded_bits_final));
    
    %decoded_bits_final.' ~= bits(1:length(decoded_bits_final));
    % Bit error rate (0 means 0% of errors, 1 means 100% of errors) 
    %bit_errors = b ~= a;
    BER(1,ii) = mean(bit_errors);
        % Finding out which bits were estimated incorrecly:
end

figure (111)
semilogy(SNR, 1/log2(M)*BER, 'LineWidth', 3); 
hold on; 
%semilogy(SNR, (1/log2(M))*P_sym_error, 'r--', 'LineWidth', 2); 
title('Bit error rate') 
xlabel('SNR [dB]') 
ylabel('BER [log scale]') 
legend('Simulated SER','Theoretical symbol error probability');

