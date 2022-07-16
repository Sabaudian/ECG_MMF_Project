function [drift] = GenDrift(Data,L)
% Generate Random values to define a random drift 
%   Argoments:
%   - Data: number of data in the dataset
%   - L: length of the signals
%   Method: r = a + (b-a).*rand(N,1) 
%           to obtain random value 
%           between the interval [a b]

% Prealloc data
A = zeros(1,Data); 
N = zeros(1,Data);
minDriftOffset = zeros(1,Data); 
maxDriftOffset = zeros(1,Data);
drift = zeros(Data,L);

% Vaiable for random genreation of data - Values ​​can be changed
% Amplitude range
aA = 0.1; bA = 0.8;
% Sinusoid period range
aN = 1; bN = 3;
% Slanted Line slope range
aMin = 0; bMin = 0.5;
aMax = 0.5; bMax = 3;
% Generate random values
for i = 1:Data
    A(i) = aA + (bA-aA).*rand(1);
    N(i) = aN + (bN-aN).*rand(1);
    minDriftOffset(i) = aMin + (bMin-aMin).*rand(1);
    maxDriftOffset(i) = aMax + (bMax-aMax).*rand(1);
    X = linspace(0,2*pi,L);
    cos_wave = A(i)*cos(X./N(i));
    slantedLine = linspace(minDriftOffset(i),maxDriftOffset(i),L);
    drift(i,:) = cos_wave + slantedLine;
end

end


