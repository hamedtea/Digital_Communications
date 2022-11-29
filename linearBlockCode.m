%LINEAR BLOCK CODER
clc;
clear all;
close all;
%%
% *Transmitter end*
in=input('Enter the (n,k)linear block code value: ');
n=in(1);
k=in(2);
p=n-k;
id=eye(k);
disp('Identity matrix')
disp(id)
temp=[ones(1,k-1) zeros(k-length(ones(1,k-1)))];
parity=[temp;];
%parity bit calculation
%for example (6,3)
%p1 = d1 xor d2
%p2 = d2 xor d3
%p3 = d3 xor d1
for i=1:n-k-1
    temp1=temp(k);
    for j=k:-1:2
        temp(j)=temp(j-1);
    end
    temp(1)=temp1;
    parity=[parity;temp]
end
parity=parity';
disp('Parity matrix')
disp(parity)
generator=[parity id];
disp('Generator matrix at transmitter end')
disp(generator)
data=0:power(2,k)-1;
data=dec2bin(data)-48;
disp('Valid data words')
disp(data)
%channel code = dataword * generator matrix
tcode=[];
for i=1:length(data)
    tcode=[tcode;mod(data(i,:)*generator,2)];
end
disp('Valid codewords ready for transmission')
disp(tcode)
%%
% *Reciever end*
rcode=input('Enter a recieved codeword: ');
check=[eye(p) parity'];
disp('Parity check matrix')
disp(check)
check=check';
disp('Transpose of Parity check matrix')
disp(check)
syn=mod(rcode*check,2);
disp('Syndrome vector in received codeoword')
disp(syn)
if syn == 0
    disp('The recieved codeword has no error in it');
    validcode=rcode;
else
    pos=0;
    for m=1:length(check)
        if(syn(1,:)==check(m,:))
            pos=m;
        end
    end
    zero=zeros(1,n);
    zero(pos)=~zero(pos);
    errorpat=zero;
    disp('Error pattern possible in recieved codeword')
    disp(errorpat);
    validcode=xor(rcode,errorpat);
    fprintf('The valid code word after detecting error at position %d and correcting it.\n',pos);
    disp(validcode)
end
%decoding data word from recieved codeword
disp('Decoded data from recieved codeword');
disp(validcode((n-k+1):n));
%minimum hamming distance
weights=[];
[r,c]=size(tcode);
for b=1:r
    weight=0;
    for o=1:length(tcode(b,:))
        if(tcode(b,o) == 1)
            weight=weight+1;
        end
    end
    weights=[weights weight];
end
[w,v]=min(weights);
weights(v)=[];
disp('Minimum Hamming Distance');
disp(min(weights));
%error detection and correction capability
disp('Error detection capability');
disp(min(weights)-1);
disp('Error correction capability');
cap=(min(weights)-1)/2;
disp(floor(cap));