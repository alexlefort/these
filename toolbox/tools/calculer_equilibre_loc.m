function NT1eq = calculer_equilibre_loc(Vs,param)

[X,U] = init_etats;
X.U = Vs;

NT1max = 100;
NT1eq = NT1max/2;
NT1min = 0.01;

for ii=1:100
    U.NT1 = NT1eq;
    [Y,DX] = equations_modele(X,U,param,0.1);
    if (DX.U < 0)
        NT1min = NT1eq;
        NT1eq = (NT1eq+NT1max)/2;
    elseif (DX.U > 0)
        NT1max = NT1eq;
        NT1eq = (NT1eq+NT1min)/2;
    end
end
