function [output] = opening(input,strel)
% Opening operation on a 1-D signal 
%   INPUT:
%   input = signal in input
%   strel = structuring element in input

    % OPENING = INPUT erosion STREL dilatation STREL
    output =dilatation(erosion(input,strel),strel);
end