M = csvread("acp_result.txt");
M = M(:,1:end-1);
[n,m] = size(M);
Mmean = M;

for ii=1:m
    Mmean(:,ii) = Mmean(:,ii) - mean(M(:,ii));
    Mmean(:,ii) = Mmean(:,ii)./sqrt(var(Mmean(:,ii)));
end

[coeff,score,latent,tsquared,explained,mu] = pca(Mmean);