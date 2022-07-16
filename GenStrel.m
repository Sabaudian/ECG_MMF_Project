function [strel] = GenStrel(N)
% Genstrel generate a structuring element
% ARGOMENT:
%   - N: define the length of strel.
%        (2*N)-1 is the dimension of the strel
 
% Gamma control the shape of strel
gamma = zeros(1,N); % Prealloc

% Height parameter, alias the peak value of strel(n) 
h = 5;
% Structuring Element
strel = zeros(1,2*N-1); % Prealloc

% generate random value for gamma [0.1-0.99]
for i = 1:length(gamma)
   gamma(i) = 0.01 + (0.99-0.01).*rand(1);
end

% Sort the first half of the array
gamma = sort(gamma);

% fill the array from index 1 to N
for n = 1:N
    strel(n) = floor(h*gamma(n));
    strel(N) = h;
end

% fill the array from index N+1 till the end
for n = N+1:2*N-1
    strel(n) = strel(2*N-n);
end

% Condition to build a correct structuring element 
if length(strel) == 1 
    strel(1) = 1;
end

for i=1:length(strel(N))
    if strel(1) > 1
       strel(1) = 0;
       strel(length(strel)) = 0;
    end
end

for i = 2:length(strel)-1
    if strel(i) < 1
        strel(i) = 1;
        strel(length(strel)-1) = 1;
    end
end

end

% %% Gen strel with only 1 element for test
% strel = zeros(1,2*N-1); % Prealloc
% % Fill the array from index 1 to N
% for n = 1:N
%     strel(n) = 1;
% end
% % Fill the array from index N+1 till the end
% for n = N+1:2*N-1
%     strel(n) = strel(2*N-n);
% end