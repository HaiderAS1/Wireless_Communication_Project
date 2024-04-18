clc; clear all;
%% 2x2, QPSK, Uncoded SM, ZF vs MMSE
Nt=2; Nr=2; % #Tx,#RX
codebook = [1+1i 1-1i -1+1i -1-1i]; %QPSK as modulation , Note Es=1+1=2;
Es=mean(abs(codebook).^2); %Symbol Energy
Eb = Es/2; %Bit Energy = Symbol Energy / log2(M);
SNR_dB=-5:2:30; No = Eb*10.^(-SNR_dB/10); %Creating Noise
BER_ZF = zeros(1,length(SNR_dB));
BER_MMSE = zeros(1,length(SNR_dB));
for i=1:length(BER_ZF)
    dist = zeros(1,length(codebook));
    n_error_ZF = 0;
    n_error_MMSE = 0;
    n_symbols = 0;
    while n_symbols <50000 %Monte Carlo Simulation Limit 5000
        s_1 = codebook(randi(4)); %1st Symbol to be transmitted
        s_2 = codebook(randi(4)); %2nd Symbol to be transmitted
        S = (1/sqrt(Nt)) * [s_1;s_2]; %Transmit Symbol
        H = (randn(Nr,Nt) + 1i*randn(Nr,Nt))/sqrt(2) ; %Channel Matrix
        n = (randn(Nr,1) + 1i*randn(Nr,1)) * sqrt(No(i)/2) ; % Noise in each time step
        y = H*S + n; %Recieved Signal Vector
        W_ZF = (inv(H'*H))*H';
        W_MMSE = (inv(H'*H+No(i)*eye(Nt,Nt)))*H';
        z_ZF = W_ZF * y; %ZF Detection
        z_MMSE = W_MMSE * y; %MMSE Detection
        for j=1:length(codebook)
            dist1_ZF(j) = abs(z_ZF(1)-codebook(j));
            dist2_ZF(j) = abs(z_ZF(2)-codebook(j));
            dist1_MMSE(j) = abs(z_MMSE(1)-codebook(j));
            dist2_MMSE(j) = abs(z_MMSE(2)-codebook(j));
        end
        [~,idx1_ZF]=min(dist1_ZF);[~,idx2_ZF]=min(dist2_ZF);
        [~,idx1_MMSE]=min(dist1_MMSE);[~,idx2_MMSE]=min(dist2_MMSE);
        s1_detected_ZF = codebook(idx1_ZF); s2_detected_ZF = codebook(idx2_ZF);
        s1_detected_MMSE = codebook(idx1_MMSE); s2_detected_MMSE = codebook(idx2_MMSE);
        if s1_detected_ZF ~= s_1
            n_error_ZF = n_error_ZF+1;
        end
        if s2_detected_ZF ~= s_2
            n_error_ZF = n_error_ZF +1;
        end
        if s1_detected_MMSE ~= s_1
            n_error_MMSE = n_error_MMSE+1;
        end
        if s2_detected_MMSE ~= s_2
            n_error_MMSE = n_error_MMSE +1;
        end
        n_symbols = n_symbols + 2;
    end
    SER_ZF = n_error_ZF / n_symbols;
    BER_ZF(i) = SER_ZF / 2; %Assuming one bit error per symbol
    SER_MMSE = n_error_MMSE / n_symbols;
    BER_MMSE(i) = SER_MMSE / 2; %Assuming one bit error per symbol
end
%%
semilogy(SNR_dB,BER_ZF,'b-d');
xlabel(' SNR (dB) ');
ylabel(' Bit Error Rate (BER) ');
hold on;
semilogy(SNR_dB,BER_MMSE,'r-s');