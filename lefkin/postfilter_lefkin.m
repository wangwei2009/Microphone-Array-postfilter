function [ z,Pxii_pre,Pxij_pre] = postfilter_lefkin( x,ybf,bf_weights,fs,N_FFT,frameLength,inc,Fvv,Pxii_pre,Pxij_pre,angle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author: wangwei
%  Data  : 6/15/2017
%  refer to:
%  "An Optimum Microphone Array Post-Filter for Speech Applications"
%  "A generalized estimation approach for linear and nonlinear
%   microphone array post-filters"
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
window = sqrt(hann(frameLength+1));
window = window(1:frameLength);
win_scale = 1;


P_len = N_FFT/2+1;
Pxii = zeros(N,P_len);
Pssnn = zeros(1,P_len);
% Pxii_pre = ones(N,P_len);
Pssnn_pre = ones(1,P_len);

Pxij = zeros((N*N-N)/2,P_len);
% Pxij_pre = ones((N*N-N)/2,P_len);
Pxij_curr = ones((N*N-N)/2,P_len);

Pss = zeros(1,P_len);
Pvv = zeros(1,P_len);
H = zeros(1,N_FFT/2+1);

Pss_e = zeros((N*N-N)/2,P_len);
Pvv_e = zeros((N*N-N)/2,P_len);

z = zeros(1,length(x(:,1)));
znoise = zeros(1,length(x(:,1)));
zspeech = zeros(1,length(x(:,1)));
t = 1;
%%
tic

f = 0:fs/256:fs/2;
w = 2*pi*fs*(0:N_FFT/2)/N_FFT;

M = N;

alpha = 0.8;

%%
Inc = inc; 
%   Wiener post-filter transfer function
%           Pss
%   h = -------------
%       Pss  +  Pnn
%
for p = 1:Inc:length(x(:,1))-N_FFT
% for p = 1:Inc:length(x(:,1))-N_FFT
    for i = 1:N
        Xi = fft(x(p:p+frameLength-1,i).*window,N_FFT);
        Pxii_curr = Xi.*conj(Xi);%abs(Xi).^2;
        % eq.11
        Pxii(i,:) = alpha*Pxii_pre(i,:)+(1-alpha)*Pxii_curr(1:N_FFT/2+1).';      
    end
    Pxii_pre = Pxii;
    Pssnn = sum(Pxii)/N;
    for i = 1:N-1
        for j = i+1:N
            
            Xi = fft(x(p:p+frameLength-1,i).*window,N_FFT).';
            Xj = fft(x(p:p+frameLength-1,j).*window,N_FFT).';
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
            Pvv_e(t,:) = (0.5*((Pxii(i,:)+Pxii(j,:)))-(real(Pxij(t,:))))...
                         ./...
                         (ones(1,P_len) - real(Fvv(:,i,j)'));
             t = t+1;
        end
    end
    Pxij_pre = Pxij;
    t = 1;
    % eq.23 
    % take the average of multichanel signal to improve robustness
    Pss = sum(Pss_e)*2/(N*N-N)*1; 
    Pvv = sum(Pvv_e)*2/(N*N-N); 
    for k = 1:N_FFT/2+1
        gamma = bf_weights(:,k)'*squeeze(Fvv(k,:,:))*bf_weights(:,k);
        H(k) = Pss(k)/(Pss(k)+1*gamma*Pvv(k));
    end
    
    % handle the indeterminite soulution when MSC¡Ö1
    % Pss(Pss<0) = 1e-3;
    
    % eq.23 
    % calculate the frequency domain filter coefficient
    beta0 = 1;
%     W_e = real(Pss)./(Pss+beta0*(Pssnn-Pss));
    W_e = real(Pss)./Pssnn;

%     W = [W_e,conj(fliplr(W_e(2:N_FFT/2)))];
    W = [H,conj(fliplr(H(2:N_FFT/2)))];
    
    % transfor the signal to frequency domain
    Xds = fft([ybf(p:p+frameLength-1)'].*window',N_FFT);
       
    % filter the signal 
    DS_filtered = W.*(Xds);
    
    % get the time domain signal
    iX = ifft(DS_filtered);
    s_est = iX(1:N_FFT);
    
    % keep the signal
    z(p:p+N_FFT-1) = z(p:p+N_FFT-1) + s_est.*window';

end

% z = z/max(abs(z));


end

