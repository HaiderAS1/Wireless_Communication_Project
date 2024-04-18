clc; clear all;
%% 2x2, QPSK, No TX CSI, No STBC
Nt=2; Nr=2; % #Tx,#RX
codebook = [1+1i 1-1i -1+1i -1-1i]; %QPSK as modulation , Note Es=1+1=2;
Es=mean(abs(codebook).^2); %Symbol Energy
Eb = Es/2; %Bit Energy = Symbol Energy / log2(M);
SNR_dB=-5:1:12; No = Eb*10.^(-SNR_dB/10); %Creating Noise
BER = zeros(1,length(SNR_dB));
for i=1:length(BER)
    dist = zeros(1,length(codebook));
    n_error = 0;
    n_symbol = 0;
    while n_symbol <50000*2 %Monte Carlo Simulation Limit 5000
        s = codebook(randi(4)); % Symbol to be transmitted
        S = (1/sqrt(Nt)) * ones(Nt,1) * s; %Same Symbol transmitted
        H = (randn(Nr,Nt) + 1i*randn(Nr,Nt))/sqrt(2) ; %Channel Matrix
        n = (randn(Nr,1) + 1i*randn(Nr,1)) * sqrt(No(i)) ; % Noise in each time step
        y = H*S + n; %Recieved Signal Vector
        heff1 = (1/sqrt(Nt))*sum(H(1,:));
        heff2 = (1/sqrt(Nt))*sum(H(2,:));
        heff_vec=[heff1;heff2];
        z = heff_vec'*y;
        for j=1:length(codebook)
            factor = sum(abs(heff_vec).^2);
            dist(j) = abs(z-factor*codebook(j)); % Distance of recieved vector with code
        end
        [~,idx]=min(dist); % Maximum Likelihood Detection
        s_detected = codebook(idx);
        if s_detected ~= s
            n_error = n_error+1;
        end
        n_symbol = n_symbol + 1;
    end
    SER = n_error / n_symbol;
    BER(i) = SER / 2; %Assuming one bit error per symbol
end
%%
semilogy(SNR_dB,BER,'r--');
xlabel(' SNR (dB) ');
ylabel(' Bit Error Rate (BER) ');
hold on;
%% 2x2, QPSK, No TX CSI, STBC - Alamouti
Nt=2; Nr=2; % #Tx,#RX
codebook = [1+1i 1-1i -1+1i -1-1i]; %QPSK as modulation , Note Es=1+1=2;
Es=mean(abs(codebook).^2); %Symbol Energy
Eb = Es/2; %Bit Energy = Symbol Energy / log2(M);
SNR_dB=-5:1:12; No = Eb*10.^(-SNR_dB/10); %Creating Noise
BER = zeros(1,length(SNR_dB));
for i=1:length(BER)
    n_error = 0;
    n_symbols = 0;
    dist1 = zeros(1,length(codebook));
    dist2 = zeros(1,length(codebook));
    while n_symbols < 50000*4
        s_1 = codebook(randi(4)); % 1st Symbol of block to be transmitted
        s_2 = codebook(randi(4)); % 2nd Symbol of block to be transmitted
        S_1 = (1/sqrt(Nt))*[s_1;s_2]; % Transmit vector of first symbol period
        S_2 = (1/sqrt(Nt))*[-conj(s_2);conj(s_1)]; % Transmit vector of second symbol period
        H = (randn(Nr,Nt) + 1i*randn(Nr,Nt))/sqrt(2) ; %Channel Matrix
        n1 =(randn(Nr,1) + 1i*randn(Nr,1)) * sqrt(No(i)) ; % Noise in first time step
        n2 =(randn(Nr,1) + 1i*randn(Nr,1)) * sqrt(No(i)) ; % Noise in second time step
        y1 = H*S_1 + n1;
        y2 = H*S_2 + n2;
        Y = [y1;conj(y2)]; %Recieved Vector
        Heff = [H; conj(H(1,2)) -conj(H(1,1)); conj(H(2,2)) -conj(H(2,1))];
        z = Heff' * Y;
        for j=1:length(codebook)
            D = Heff'*Heff;
            dist1(j) = abs(z(1)-D(1,1)*codebook(j)); %Take Eucledian distance
            dist2(j) = abs(z(2)-D(2,2)*codebook(j));
        end
        [~,idx1] = min(dist1); %Find the minimum index
        [~,idx2] = min(dist2);
        s1_detected = codebook(idx1); s2_detected = codebook(idx2);
        if s1_detected ~= s_1
            n_error = n_error+1;
        end
        if s2_detected ~= s_2
            n_error = n_error+1;
        end
        n_symbols = n_symbols + 2;
    end
    SER = n_error / n_symbols;
    BER(i) = SER / 2; %Assuming one bit error per symbol
end
%%
semilogy(SNR_dB,BER,'b--');
xlabel(' SNR (dB) ');
ylabel(' Bit Error Rate (BER) ');
hold on;
%% 2x1, QPSK, TX CSI, Dominant EigenMode Transmission
Nt=2; Nr=2; % #Tx,#RX
codebook = [1+1i 1-1i -1+1i -1-1i]; %QPSK as modulation , Note Es=1+1=2;
Es=mean(abs(codebook).^2); %Symbol Energy
Eb = Es/2; %Bit Energy = Symbol Energy / log2(M);
SNR_dB=-5:1:12; No = Eb*10.^(-SNR_dB/10); %Creating Noise
BER = zeros(1,length(SNR_dB));
for i=1:length(BER)
    dist = zeros(1,length(codebook));
    n_error = 0;
    n_symbol = 0;
    while n_symbol <50000*4 %Monte Carlo Simulation Limit 5000
        s = codebook(randi(4)); % Symbol to be transmitted
        H = (randn(Nr,Nt) + 1i*randn(Nr,Nt))/sqrt(2) ; %Channel Matrix
        [U,D,V] = svd(H); %SVD of Channel Matrix
        w = sqrt(Nt)*V(:,1); %Precode Vector
        S = w*s/sqrt(Nt); %Transmitted Vector
        n = (randn(Nr,1) + 1i*randn(Nr,1)) * sqrt(No(i)) ; % Noise in each time step
        y = H*S + n; %Recieved Signal Vector
        g = U(:,1); %Postcode Vector
        z = g'*y;
        for j=1:length(codebook)
            heff = D(1,1);
            dist(j) = abs(z-heff*codebook(j)); % Distance of recieved vector with code
        end
        [~,idx]=min(dist); % Maximum Likelihood Detection
        s_detected = codebook(idx);
        if s_detected ~= s
            n_error = n_error+1;
        end
        n_symbol = n_symbol + 1;
    end
    SER = n_error / n_symbol;
    BER(i) = SER / 2; %Assuming one bit error per symbol
end
%%
semilogy(SNR_dB,BER,'k--');
xlabel(' SNR (dB) ');
ylabel(' Bit Error Rate (BER) ');
hold on;