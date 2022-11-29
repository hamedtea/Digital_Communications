%channel estimation by using training sequences  
N_ref_symbols = 30; %number of reference symbols used for channel estimation 
estimate_length = 7; %this defines how long is the channel estimate's impulse response 
A_conv=convmtx(st(2:N_ref_symbols+1).',estimate_length); % Convolution matrix 
p_LS=((A_conv'*A_conv)\(A_conv'))*qk(1:size(A_conv,1)).'; % LS solution
[H_p_LS,f_p_LS]=freqz(p_LS,1,-Rs/2:Rs/400:Rs/2,Rs);

%plotting estimated channel amplitude response + actual amplitude response of the channel
figure(8); 
hold on  
plot(f_ch/1e6,20*log10(abs(H_ch)),'b')
plot(f_p_LS/1e6,20*log10(abs(H_p_LS)),'r')
legend('Channel','LS Channel Estimatation');
figure(9); hold on;  
stem(-L1:L2,abs(h),'b');
stem(-L1:length(p_LS)-L1-1,abs(p_LS),'r');
legend('Channel','LS channel estimate'); 
title('Absolute values of the impulse responses'

M = [1 2 3; 4 5 6; 7 8 9]
len = 10;                     % Number of rows in new matrix
V = randi(size(M,1), len, 1) % Random row indices
N = M(V,:)



b=rand(40,1).';
for i=1:7:20
    j=1;
    a=b(i:1:i+4);
end





        if s==[0 0 0]
            e=[0 0 0 0 0 0 0];
            y=bitxor(r,e);
        elseif s==[1 0 1]
            e=[0 0 0 0 0 0 1];
            y=bitxor(r,e);
        elseif s==[0 1 1]
            e=[0 0 0 0 0 1 0];
            y=bitxor(r,e);
        elseif s==[1 1 0]
            e=[0 0 0 0 1 0 0];
            y=bitxor(r,e);
        elseif s==[0 0 1]
            e=[0 0 0 1 0 0 0];
            y=bitxor(r,e);
        elseif s==[0 1 0]
            e=[0 0 1 0 0 0 0];
            y=bitxor(r,e);
        elseif s==[1 0 0]  
            e=[0 1 0 0 0 0 0];
            y=bitxor(r,e);
        elseif s==[1 1 1]
           e=[1 1 1 1 1 1 1];
           y=bitxor(r,e);
        else 
            y=r;
        end
        decoded_bits=[decoded_bits(ii,:) y(1:k)];
    end
end
end



N   = 10;
mat1 = rand(N);
vec1 = zeros(N*N,1);
for i=1:N
    for j=1:N
        vec1((i-1)*N + j) = mat1(i,j);
    end
end



%[numofrows,numofcolumns] = size(encoded_bits);
%encoded_bits_final = zeros(numofrows*numofcolumns,1);
%for i=1:numofrows
%    for j=1:numofcolumns
%        encoded_bits_final(i) = encoded_bits(i,j);
%    end
%end