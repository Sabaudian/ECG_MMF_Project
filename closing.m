function [output] = closing(input,strel)
% Closing operation on a 1-D signal 
%   INPUT:
%   input = signal in input
%   strel = structuring element in input
    
    % OPENING = INPUT dilatation STREL erosion STREL
    output = erosion(dilatation(input,strel),strel);
end