function [output] = dilatation(input,strel)
% Dilatation operation on a 1-D signal
%   INPUT:
%   insputSignal = the signal
%   strel = structuring element in input

    % Get sizes
    M = length(strel);
    N = length(input);
    
    % Strel being odd
    if mod(M,2) == 0
        error('strel lenght must be odd')
    end

    % Pad signal
    hw = floor(M/2);
    input = padarray(input,[0 hw],'replicate','both');
    
    % Perform dilatation
    output = zeros(1,N);
    for n = (1:N)+hw % hw+1:hw+N
         output(n-hw) = max(input(n-hw:n+hw) + strel);
    end
end