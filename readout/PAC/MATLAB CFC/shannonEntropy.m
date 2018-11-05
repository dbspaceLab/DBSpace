function entropy = shannonEntropy(inputVector)
%SHANNONENTROPY Determines statistical entropy of input probability vector
%   USAGE: entropy = shannonEntropy(inputVector)
%   Calculates SUM of all P(i) * log ( P(i) ) and outputs average
%   probability of all elements (P) of inputVector
    entropy = -sum(inputVector.*log(inputVector));
end

