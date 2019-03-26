function [ Fvv ] = GenNoiseMSC(M,N,fs,r)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明

N_FFT = N;
c = 340;
f = 0:fs/N_FFT:fs/2;
Fvv = zeros(N_FFT/2+1,M,M);
k_optimal = 1;
for i = 1:M
    for j = 1:M   
        if i == j
            Fvv(:,i,j) = ones(N_FFT/2+1,1);
        else
            dij = abs(i-j)*r;
            Fvv(:,i,j) = sin(2*pi*f*dij*k_optimal/c)./(2*pi*f*dij*k_optimal/c);Fvv(1,i,j) = 0.998;%T(2) = 0.996;
        end
        index = find(Fvv(:,i,j)>0.7);
        if(size(index,1)>0)
            Fvv(index,i,j)=0.7;
        end
    end
end

end

