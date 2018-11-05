function MI = ModulationIndex(P)
%MODULATION INDEX determines the MI for an input Probability Vector P
%   USAGE: MI = ModulationIndex(P)
%   Divides the Kullback-Liebler Divergence of P with respect to the
%   uniform distribution by natural log of the length of P which bounds the
%   output between 0 and 1
    MI = KLDiv(P)/length(P);
end

