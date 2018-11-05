function KullbackLieblerDivergence = KLDiv(varargin)
%KLDIV Calculates the Kullback-Liebler Divergence of the probability vector
%P from the probability vector Q. 
%   USAGE KullbackLieblerDivergence = KLDiv(P,Q)
%   If Q is not given, it is assumed to be the uniform distribution of the 
%   length of P. This function accepts two inputs methods: KLDiv(P,Q) or
%   KLDiv(P) with Q implied as aforementioned.
switch nargin
    case 1
        P = varargin{1}; entropyP = shannonEntropy(P);
        KullbackLieblerDivergence = log(length(P))-entropyP;
    case 2
        P = varargin{1}; Q = varargin{2};
        if length(P) ~= length(Q)
            error('Lengths of P and Q probability vectors must be equal');
        end
        KullbackLieblerDivergence = sum(P.*log(P./Q));
    case 0
        error('Too few inputs: Needs atleast one input vector');
    otherwise
        error('Unexpected Inputs: Too many inputs');
end

