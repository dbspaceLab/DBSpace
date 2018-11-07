function zScore = zScoredMVL(MultChannIn)
%ZSCOREDMVL Give a z-score to the mean vector based on PCA
%   USAGE: zScore = zScoredMVL(MultChannIn)
%   Each row of MultChannIn is a separate channel, each column represents a
%   synchronous sampling of all channels at a time point.  
%   zScore is the z-score of the mean column vector in MultChannIn
    [XPrinVals, XPrinComps] = PrinCompAnal(MultChannIn);
    meanVect = [mean(MultChannIn(1,:)); mean(MultChannIn(2,:))];
    theta = acos(meanVect'*XPrinComps(:,1)/norm(meanVect));
    R = sqrt((sqrt(XPrinVals(1))*cos(theta))^2+(sqrt(XPrinVals(2))*sin(theta))^2);
    zScore = norm(meanVect)/R;
end

