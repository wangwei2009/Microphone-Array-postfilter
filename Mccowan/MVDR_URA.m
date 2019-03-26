function [ MVDRout,H,DI,WNG] = MVDR_URA( x,fs,N,frameLength,inc,r,angle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency domain MVDR beamformer
%   
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
if N_argout > 2
   disp('Calculating isotropic noise MSC')
   Fvv = GenNoiseMSC(size(x,2),N,fs,r);
end
c = 340;
Nele = size(x,2);
omega = zeros(frameLength,1);
H = ones(N/2+1,Nele)';

theta = angle(2);%90*pi/180; %固定一个俯仰角
gamma = [0 90 180 270]'*pi/180;%麦克风位置
% gamma = [30 90 150 210 270 330]*pi/180;%麦克风位置
tao = -1*r*cos(theta)*cos(angle(1)-gamma)/c;     %方位角 0 < angle <360
yds = zeros(length(x(:,1)),1);
x1 = zeros(size(x));

DI = zeros(1,N/2+1);
WNG = zeros(1,N/2+1);

M = Nele;

alpha = 0.9;
Rvv = zeros(M,M,N/2+1);
for i = 1:inc:length(x(:,1))-frameLength


    Z = fft(x(i:i+frameLength-1,:).*hamming(frameLength)).';
    
    for k = 1:N/2+1
        omega(k) = 2*pi*(k-1)*fs/N;    
        
        % propagation vector
        a = exp(-1j*omega(k)*tao);
        
        if(i<10000)
            Rvv(:,:,k) = squeeze(alpha*Rvv(:,:,k))+(1-alpha)*Z(:,k)*Z(:,k)';
        end      
        Fvv_k = (squeeze(Rvv(:,:,k))+1e-1*eye(Nele));% Diagonal loading
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
    x_fft = conj(H).*Z(:,1:N/2+1);
    yf = sum(x_fft);
    Cf = [yf,conj(fliplr(yf(2:N/2)))];
    
    % 恢复延时累加的信号
    yds(i:i+frameLength-1) = yds(i:i+frameLength-1)+(ifft(Cf)');
    
end
MVDRout = real(yds);  
if(N_argout>=3)
   DI = pow2db(real(DI));
end
if(N_argout>=4)
    WNG = pow2db(real(WNG));
end

end

