function [Xeq,Ueq] = calculer_equilibre(Vs, param, r)

deg = pi/180;
%% Initialisation des conditions d'equilibre

[Xeq,Ueq] = init_etats;
Xeq.U = Vs;

%% Révision du point d'equilibre en immersion

if (~isfield('b20',r))
    r.b10 = 1e5;
    r.b20 = 1e-5;
    r.epsilon0 = 1;
    r.theta0 = 1;
end

[~,PenteEq,Xeq.Theta,Ueq.Beta1,Ueq.Beta2] = EquilibreLongitudinal(0.0, 20.0*deg, Vs, param, r, 0);
Xeq.W = Vs*tan(Xeq.Theta - PenteEq);

%% Calcul du point d'equilibre en vitesse (NT1Eq pour U)

Ueq.NT1 = calculer_equilibre_loc(Vs, param);
