function [noise] = GenNoise(N,L)
%Generate randon Gaussian noise
%   INPUT:
%   - N: number of data in the dataset (row)
%   - L: length of the signals (column)
%   Method: r = a + (b-a).*rand(N,1) 
%           to obtain random value 
%           between the interval [a b]

% Prealloc data
E = zeros(1,N); % 
Sigma1 = zeros(1,N); % Standard deviation
Sigma2 = zeros(1,N); % Standard deviation
snr = zeros(1,N); % SNR
noise = zeros(N,L); % output noise
% Epsilon Range
aE = 0.1; bE = 0.5;
% As Sigma1 and Sigma2 increases, the noise amplitude increases
% Standard Deviation Range
aS1 = 0.01; bS1 = 0.03;
% SNR range
aSnr = 0.5; bSnr = 1.5;
% Generate random variable to gen. noise
for i = 1:N
    E(i) = aE + (bE-aE).*rand(1);
    Sigma1(i) = aS1 + (bS1-aS1).*rand(1);
    aS2 = 2*Sigma1(i); bS2 = 20*Sigma1(i);
    Sigma2(i) = aS2 + (bS2-aS2).*rand(1);
    snr(i) = aSnr + (bSnr - aSnr).*rand(1);
end
% Generate noise
for i = 1:N
    noise(i,:) = snr(i)*((1-E(i))*Sigma1(i).*randn(L,1) + E(i)*Sigma2(i).*randn(L,1));
end
% Sort the noise ascending order
noise = sort(noise); 

end