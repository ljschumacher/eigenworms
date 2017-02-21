function D_KL = kldiv(P,Q)
% KLDIV calculates Kullblack-Leibler divergence between two n-dimensionsl
% probability distributions P and Q

assert(all(size(P)==size(Q)))

D_KL = sum(P.*log(P./Q),'omitnan');
while length(D_KL)>1
    D_KL = sum(D_KL);
end
end