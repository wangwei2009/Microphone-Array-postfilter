function [ z,Pxii_pre,Pxij_pre] = post_filter_func( x,DS,fs,N_FFT,frameLength,inc,Fvv,Pxii_pre,Pxij_pre,angle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  General post-filter based on noise filed coherence 
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

N = size(x,2);        %Channels
M = N;
% angle = [90,0]/180*pi;
%% 


% periodic window for spectral analysis
window = hamming(frameLength+1);
win = window(1:frameLength);
win_scale = 1.08;


P_len = N_FFT/2+1;
Pxii = zeros(N,P_len);
Pssnn = zeros(1,P_len);
% Pxii_pre = ones(N,P_len);
Pssnn_pre = ones(1,P_len);

Pxij = zeros((N*N-N)/2,P_len);
% Pxij_pre = ones((N*N-N)/2,P_len);
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

alpha = 0.5;

%%
Inc = 128; 
%   Wiener post-filter transfer function
%           Pss
%   h = -------------
%       Pss  +  Pnn
%
for p = 1:Inc:length(x(:,1))-N_FFT
% for p = 1:Inc:length(x(:,1))-N_FFT
    for i = 1:N
        Xi = fft(x(p:p+frameLength-1,i).*win,N_FFT);
        Pxii_curr = Xi.*conj(Xi);%abs(Xi).^2;
        % eq.11
        Pxii(i,:) = alpha*Pxii_pre(i,:)+(1-alpha)*Pxii_curr(1:N_FFT/2+1).';      
    end
    Pxii_pre = Pxii;
    Pssnn = sum(Pxii)/N;
    for i = 1:N-1
        for j = i+1:N
            
            Xi = fft(x(p:p+frameLength-1,i).*win,N_FFT).';
            Xj = fft(x(p:p+frameLength-1,j).*win,N_FFT).';
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
    beta0 = 1;
%     W_e = real(Pss)./(Pss+beta0*(Pssnn-Pss));
    W_e = real(Pss)./Pssnn;

    W = [W_e,conj(fliplr(W_e(2:N_FFT/2)))];
    
    % transfor the signal to frequency domain
    Xds = fft([DS(p:p+frameLength-1)'].*win',N_FFT);
       
    % filter the signal 
    DS_filtered = W.*(Xds);
    
    % get the time domain signal
    iX = ifft(DS_filtered);
    s_est = iX(1:N_FFT);
    
    % keep the signal
    z(p:p+N_FFT-1) = z(p:p+N_FFT-1) + s_est;

end

z = z/win_scale;




end

