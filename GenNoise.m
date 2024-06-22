function [noise] = GenNoise(N,L)
% Generate random Gaussian noise as a mixture model
%   INPUT:
%   - N: number of data in the dataset (row)
%   - L: length of the signals (column)
%   Method: r = a + (b-a).*rand(N,1) 
%           to obtain random value 
%           between the interval [a b]

% Prealloc data
E = zeros(1,N); % Epsilon (weight)
Sigma1 = zeros(1,N); % Standard deviation of G1
Sigma2 = zeros(1,N); % Standard deviation of G2
noise = zeros(N,L); % output noise

% Epsilon Range
aE = 0.1; bE = 0.5;

% As Sigma1 and Sigma2 increase, the noise amplitude increases
% Standard Deviation Range for Sigma1
aS1 = 0.01; bS1 = 0.03;

% Generate random variables to gen. noise
for i = 1:N
    E(i) = aE + (bE-aE).*rand(1); % Generate epsilon
    Sigma1(i) = aS1 + (bS1-aS1).*rand(1); % Generate Sigma1
    
    % Standard Deviation Range for Sigma2
    aS2 = 2*Sigma1(i); bS2 = 20*Sigma1(i);
    Sigma2(i) = aS2 + (bS2-aS2).*rand(1); % Generate Sigma2
end

% Generate noise
for i = 1:N
    % Generate Gaussian noise with standard deviation Sigma1
    G1 = Sigma1(i) .* randn(L,1);
    % Generate Gaussian noise with standard deviation Sigma2
    G2 = Sigma2(i) .* randn(L,1);
    % Combine them according to the mixture model
    noise(i,:) = (1 - E(i)) * G1' + E(i) * G2';
end

% Sort the noise in ascending order
noise = sort(noise, 2); % Sort along each row

end
