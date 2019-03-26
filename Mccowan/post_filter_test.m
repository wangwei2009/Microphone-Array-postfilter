%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  General post-filter based on noise filed coherence 
%    test with office noisy recordings
%    array type:4 mic URA,0.032cm
%  dependencies:
%    RIR-Generator
%    Signal-Generator
%
%  Author: wangwei
%  Data  : 6/15/2017
%  refer to:
%  "Microphone Array Post-Filter Based on Noise Field Coherence" IEEE 2003
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
% clear all;
c = 340; % speed of sound

%%
%% load recorded office noise audio
noisepath = '../sound/rec1/';
[noise1,fs] = audioread([noisepath,'Òô¹ì-2.wav']);
noise2 = audioread([noisepath,'Òô¹ì-3.wav']);
noise3 = audioread([noisepath,'Òô¹ì-4.wav']);
noise4 = audioread([noisepath,'Òô¹ì-5.wav']);
%gsc_out = audioread('output/out_gsc1.wav');
x = [noise1,noise2,noise3,noise4];

bhi = fir1(512,0.01,'high');
x = filter(bhi,1,x);
%x = x(1:length(gsc_out),:);
interval = 1:120000;
% x = x(interval,:);
N = size(x,2);        %Channels
M = N;
angle = [195,0]/180*pi;
r = 0.032;
%% 
N_FFT = 256;

window = hamming(N_FFT);


P_len = N_FFT/2+1;
Pxii = zeros(N,P_len);
Pssnn = zeros(1,P_len);
Pxii_pre = ones(N,P_len);
Pssnn_pre = ones(1,P_len);

Pxij = zeros((N*N-N)/2,P_len);
Pxij_pre = ones((N*N-N)/2,P_len);
Pxij_curr = ones((N*N-N)/2,P_len);

Pss = zeros(1,P_len);

Pss_e = zeros((N*N-N)/2,P_len);

z = zeros(1,length(x(:,1)));
znoise = zeros(1,length(x(:,1)));
zspeech = zeros(1,length(x(:,1)));
t = 1;
%%
tic

f = 0:fs/256:fs/2;
w = 2*pi*fs*(0:N_FFT/2)/N_FFT;

M = N;

alpha = 0.9;
Fvv = zeros(N_FFT/2+1,M,M);
%% Frequency domain delay-sum,time alignment
[ DelaySumOut, x1] = DelaySumURA(x,fs,N_FFT,N_FFT,N_FFT/2,r,angle);


%% estimate noise coherence function
noise = x1(1:40000,:);
% x = x(40000:end,:);
for i = 1:size(noise,2)
    for j = 1:size(noise,2)
        
        [sc,F]=mycohere(noise(:,i),noise(:,j),256,fs,hanning(256),0.75*256);
        Fvv(:,i,j) = real(sc);

        index = find(Fvv(:,i,j)>0.9);
        if(size(index,1)>0)
            Fvv(index,i,j)=0.9;
        end
%         Fvv(1,i,j) = 0.90;
    end
end
Fvv(1,:,:) = 0.90;

%% fixed wideband MVDR using pre-defined noise cross-spectral matrix
% [ MVDR_out,H,DI,WNG] = superdirectiveMVDR(x,fs,N_FFT,N_FFT,N_FFT/2,r,angle);
[ MVDR_out] = MVDR_URA(x,fs,N_FFT,N_FFT,N_FFT/2,r,angle);

x = x1;
DS = MVDR_out;
% [sc,F]=mycohere(x1(:,1),x1(:,2),256,fs,hanning(256),0.75*256);
% Fvv2 = Fvv;

%%
Inc = 128; 
k_optimal = 1;
hwt = waitbar(0,'general poster filter');
%   Wiener post-filter transfer function
%           Pss
%   h = -------------
%       Pss  +  Pnn
%
for p = 1:Inc:length(x(:,1))-N_FFT
    for i = 1:N
        Xi = fft(x(p:p+N_FFT-1,i).*window);
        Pxii_curr = Xi.*conj(Xi);%abs(Xi).^2;
        % eq.11
        Pxii(i,:) = alpha*Pxii_pre(i,:)+(1-alpha)*Pxii_curr(1:N_FFT/2+1).';      
    end
    Pxii_pre = Pxii;
    Pssnn = sum(Pxii)/N;
    for i = 1:N-1
        for j = i+1:N
            
            Xi = fft(x(p:p+N_FFT-1,i).*window).';
            Xj = fft(x(p:p+N_FFT-1,j).*window).';
            % cross-spectral
            Pxij_temp = Xi.*conj(Xj);
            % half bin
            Pxij_curr(t,:) = Pxij_temp(1:N_FFT/2+1);
            % average
            Pxij(t,:) = alpha*Pxij_pre(t,:)+(1-alpha)*Pxij_curr(t,:);
                 
            % eq.22 estimate source signal's PSD
            Pss_e(t,:) = (real(Pxij(t,:)) - 0.5*real(Fvv(:,i,j)').*(Pxii(i,:)+Pxii(j,:)))...
                         ./...
                         (ones(1,P_len) - real(Fvv(:,i,j)'));
             t = t+1;
        end
    end
    Pxij_pre = Pxij;
    t = 1;
    % eq.23 
    % take the average of multichanel signal to improve robustness
    Pss = sum(Pss_e)*2/(N*N-N); 
    
    % handle the indeterminite soulution when MSC¡Ö1
    % Pss(Pss<0) = 1e-3;
    
    % eq.23 
    % calculate the frequency domain filter coefficient
    W_e = real(Pss)./Pssnn;

    W = [W_e,conj(fliplr(W_e(2:128)))];
    
    % transfor the signal to frequency domain
    Xds = fft([DS(p:p+N_FFT-1)'].*window');
       
    % filter the signal 
    DS_filtered = W.*(Xds);
    
    % get the time domain signal
    iX = ifft(DS_filtered);
    s_est = iX(1:N_FFT);
    
    % keep the signal
    z(p:p+N_FFT-1) = z(p:p+N_FFT-1) + s_est;

    waitbar(p/(length(x(:,1))));
end
close(hwt);
toc
% z = z/max(abs(z));


