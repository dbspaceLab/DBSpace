function [PrinVals PrinComps] = PrinCompAnal(MultChannIn)
%PRINCOMPANAL Outputs Principal Component Analysis of Multichannel Input                                                                       
%   USAGE: [PrinVals PrinComps] = PrinCompAnal(MultChannIn)
%   Each row of MultChannIn is a separate channel, each column represents a
%   synchronous sampling of all channels at a time point.   
%   PrinVals is a column vector containing the eigenvalues of the
%   covariance matrix for MultChannIn.
%   PrinComps is a matrix whose columns are the principal components of the
%   MultChannIn data.
    MultChannInCov = cov(MultChannIn');
    [PrinComps, PrinValsDiag] = eig(MultChannInCov);
    PrinVals = zeros(length(PrinValsDiag),1);
    for kk = 1:length(PrinVals)
        PrinVals(kk) = PrinValsDiag(kk,kk);
    end
end

