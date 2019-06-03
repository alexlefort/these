load('transferts.mat');

M = transferts.mat_stab;

f = fields(model_sm.p0);
for ii=1:length(f)
    M = subs(M,f{ii},eval(model_sm.p0.(f{ii})));
end

rng(1234);

kpsi_min =   0;
kpsi_max =  100;
kr_min   =  -1;
kr_max   =  100;
ir_min   =   0;
ir_max   =  100;

n = 10000;
res = zeros(n,4);

syms Mev kpsi kr ir;
Mev(kpsi,kr,ir) = M;


for ii=1:n    
    kpsie = (kpsi_max - kpsi_min)*rand(1) + kpsi_min;
    kre   = (ir_max   - kr_min)*rand(1)   + kr_min;
    ire   = (ir_max   - ir_min)*rand(1)   + ir_min;
    Maux = eval(Mev(kpsie,kre,ire));
    res(ii,1) = kpsie;
    res(ii,2) = kre  ;
    res(ii,3) = ire  ;
    res(ii,4) = (max(real(eig(Maux)))<0);
    if (mod(ii,100) == 0)
        disp(ii);
        ratio = sum(res(:,4))/ii;
        disp(ratio);
    end
end

ratio = sum(res(:,4))/n;
disp(ratio);


%% Get all stabilizing controllers


jj=1;
stab = zeros(1,3);

for ii=1:n
    if res(ii,4) == 1
        stab(jj,1:3) = res(ii,1:3);
        jj = jj+1;
    end
end

