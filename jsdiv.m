function D_JS = jsdiv(P,Q)
% JSDIV calculates Jensen-Shannon divergence between two n-dimensionsl
% probability distributions P and Q

assert(all(size(P)==size(Q)))

M = (P + Q)/2;
D_JS = kldiv(P,M)/2 + kldiv(Q,M)/2;

end

