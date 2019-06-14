%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test lefkim postfilter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
% clear all;
c = 340; % speed of sound

%%
%% load recorded office noise audio
noisepath = '../an101-mtms-arrA/an101-mtms-arrA ';
% noisepath = '../sound/backward/2/';
[noise1,fs] = audioread([noisepath,'1.wav']);
noise2 = audioread([noisepath,'2.wav']);
noise3 = audioread([noisepath,'3.wav']);
noise4 = audioread([noisepath,'4.wav']);
noise5 = audioread([noisepath,'5.wav']);
noise6 = audioread([noisepath,'6.wav']);
noise7 = audioread([noisepath,'7.wav']);
noise8 = audioread([noisepath,'8.wav']);
x = [noise1,noise2,noise3,noise4,noise5,noise6,noise7,noise8];
% x = [noise2,noise3,noise4,noise5,noise6,noise7];
%% high-pass
bhi = fir1(512,0.01,'high');
x = filter(bhi,1,x);

interval = 1:120000;
% x = x(interval,:);
N = size(x,2);        %Channels
M = N;
angle = [197,90]/180*pi;
r = 0.07;
frameLength = 256;
overlap = 128;
inc = frameLength - overlap;
N_FFT = 256;

P_len = N_FFT/2+1;
Pxii_pre = ones(N,P_len);
Pxij_pre = ones((N*N-N)/2,P_len);

%% Frequency domain delay-sum,time alignment
% [ DelaySumOut, x1,DS_weight] = DelaySumURA(x,fs,N_FFT,frameLength,inc,r,angle);
Fvv_th = GenNoiseMSC(M,N_FFT,fs,r);
%% fixed wideband MVDR using pre-defined noise cross-spectral matrix
angle = [197,90]/180*pi;
% [ MVDR_out,H,DI,WNG] = superdirectiveMVDR(x1,fs,N_FFT,frameLength,inc,r,angle);
MVDRweights = ones(M,N_FFT/2+1)*0.25;
%%
[ z] = postfilter_lefkin( x,sum(x,2)/8,MVDRweights,fs,N_FFT,frameLength,inc,Fvv_th,Pxii_pre,Pxij_pre,angle);

