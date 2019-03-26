function [ MVDRout,H,DI,WNG] = superdirectiveMVDR( x,fs,N_FFT,frameLength,inc,r,angle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency domain MVDR beamformer
% [ MVDRout,H,DI,WNG] = superdirectiveMVDR( x,fs,N,frameLength,inc,r,angle)  
%      input :
%          x : input signal ,samples * channel
%          fs: sample rate
%          N : fft length,frequency bin number
%frameLength : frame length,usually same as N
%        inc : step increment
%          d : array element spacing
%      angle : incident angle
%
%     output :
%    MVDRout : MVDR output
%         x1 : presteered signal,same size as x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_argout = nargout;
DS_weight = GetDSweight(size(x,2),N_FFT,fs,r,angle);
% if N_argout > 2
   disp('Calculating isotropic noise MSC')
   Fvv = GenNoiseMSC_shift(size(x,2),N_FFT,fs,r,DS_weight);
% end

c = 340;
Nele = size(x,2);
omega = zeros(frameLength,1);
H = ones(N_FFT/2+1,Nele)';

theta = angle(2);%90*pi/180; %固定一个俯仰角
gamma = [0 90 180 270]'*pi/180;%麦克风位置
% gamma = [30 90 150 210 270 330]*pi/180;%麦克风位置
tao = -1*r*cos(theta)*cos(angle(1)-gamma)/c;     %方位角 0 < angle <360
yds = zeros(length(x(:,1)),1);
x1 = zeros(size(x));

DI = zeros(1,N_FFT/2+1);
WNG = zeros(1,N_FFT/2+1);

M = Nele;

window = hamming(frameLength+1);
window = window(1:frameLength); 
win_scale = 1.08;

for i = 1:inc:length(x(:,1))-N_FFT


    Z = fft(x(i:i+frameLength-1,:).*window,N_FFT).';
    
    for k = 1:N_FFT/2+1
        omega(k) = 2*pi*(k-1)*fs/N_FFT;    
        
        % propagation vector
        a = exp(-1j*omega(k)*tao);
            
        Fvv_k = (squeeze(Fvv(k,:,:))+1e-1*eye(Nele));% Diagonal loading

        H(:,k) =    Fvv_k\a ...
                 ./(a'/Fvv_k*a);                     % MVDR weights
             
        if N_argout>=3
        % directivity index
        DI(k) = (abs(H(:,k)'*a))^2 ...
                /(H(:,k)'*squeeze(Fvv(k,:,:))*H(:,k));
        end
        if N_argout==4
            % White Noise Gain
            WNG(k) = (abs(H(:,k)'*a))^2 ...
                    /(H(:,k)'*H(:,k));
        end
    
    end
    x_fft = conj(H).*Z(:,1:N_FFT/2+1);
    yf = sum(x_fft);
    Cf = [yf,conj(fliplr(yf(2:N_FFT/2)))];
    
    % 恢复延时累加的信号
    yds(i:i+N_FFT-1) = yds(i:i+N_FFT-1)+(ifft(Cf)');
    
end
MVDRout = real(yds)/win_scale;  
if(N_argout>=3)
   DI = pow2db(real(DI));
end
if(N_argout>=4)
    WNG = pow2db(real(WNG));
end
end

