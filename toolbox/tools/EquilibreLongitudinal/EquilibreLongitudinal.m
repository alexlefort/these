function [valid,PenteEq,ThetaEq,Beta1Eq,Beta2Eq] = EquilibreLongitudinal(PenteCom,ThetaCible,VsEst,param, reglages,equivalentX)
% calcul de la pente de trajectoire, de l'assiette et du braquage des
% barres � l'�quilibre, au plus proche de la pente de trajectoire commandee

if (~isfield('b20',reglages))
    reglages.b10 = 1e5;
    reglages.b20 = 1e-5;
    reglages.epsilon0 = 1;
    reglages.theta0 = 1;
end


param.b10      = reglages.b10       ;
param.b20      = reglages.b20       ;
param.theta0   = reglages.theta0    ;
param.epsilon0 = reglages.epsilon0  ;
param.garde    = 0.0;

% mise en forme de l'equation d'equilibre :
%  Beta1Eq = c2.ThetaEq� + c1*ThetaEq + cpi*PenteEq + c0
%  Beta2Eq = d2.ThetaEq� + d1*ThetaEq + dpi*PenteEq + d0

if (equivalentX == 1)
    param.CzB1 = (param.CzAL1 + param.CzAL2 + param.CzAL3 + param.CzAL4);
    param.CmB1 = (param.CmAL1 + param.CmAL2 + param.CmAL3 + param.CmAL4);
end

MuEst = param.mu;
XgEst = param.xg;

delta = param.CzB1*param.CmB2-param.CmB1*param.CzB2 ;
m11 = param.CmB2/delta ;
m12 = -param.CzB2/delta ;
m21 = -param.CmB1/delta ;
m22 = param.CzB1/delta ;

Vs2 = max(VsEst*VsEst,0.001) ;

k1 = 2*(MuEst-1)*param.g*param.CRemp*param.L/Vs2 ;
k2 = 2*MuEst*param.g*param.CRemp/Vs2 ;

a2 = k1/2 ;
a1 = - param.CzW ;
api = param.CzW ;
a0t = -param.Cz0 - k1 ;

b2 = -k2*XgEst/2 ;
b1 = -param.CmW + k2*param.Zg;
bpi = param.CmW ;
b0t = - param.Cm0 + k2*XgEst;

data.c2 =  m11*a2  + m12*b2 ;
data.c1 =  m11*a1  + m12*b1 ;
data.cpi = m11*api + m12*bpi ;
data.c0t =  m11*a0t  + m12*b0t ;

data.d2 =  m21*a2  + m22*b2 ;
data.d1 =  m21*a1  + m22*b1 ;
data.dpi = m21*api + m22*bpi ;
data.d0t =  m21*a0t  + m22*b0t ;


% Recherche de points d'intersection avec les bords du domaine pour savoir
% si la pente command�e est admissible
nb_racines_admissibles = 0 ; % entre 0 et 8
racines_admissibles = [0,0,0,0,0,0,0,0] ;

nb_intervalles = 1 ; % au max 5
bornes_inf_intervalles = [0,0,0,0,0] ;
bornes_sup_intervalles =  [0,0,0,0,0] ;

% Recherche des points d'intersection de l'espace solution avec
% Beta1=Beta1min
[nbsols,racines] = RacinesReellesPolynomeDegre3(0,data.c2,data.c1,(data.c0t + data.cpi*PenteCom-param.Beta1Min)) ;
for ii=1:nbsols
    beta2 = data.d2*racines(ii)*racines(ii) + data.d1*racines(ii) + data.d0t + data.dpi*PenteCom ;
    if ((abs(racines(ii)) < ThetaCible)&&(beta2 < param.Beta2Max)&&(beta2 > param.Beta2Min))
        nb_racines_admissibles = nb_racines_admissibles + 1 ;
        racines_admissibles(nb_racines_admissibles) = racines(ii) ;
    end
end

% Recherche des points d'intersection de l'espace solution avec
% Beta1=Beta1max
[nbsols,racines] = RacinesReellesPolynomeDegre3(0,data.c2,data.c1,(data.c0t + data.cpi*PenteCom-param.Beta1Max)) ;
for ii=1:nbsols
    beta2 = data.d2*racines(ii)*racines(ii) + data.d1*racines(ii) + data.d0t + data.dpi*PenteCom ;
    if ((abs(racines(ii)) < ThetaCible)&&(beta2 < param.Beta2Max)&&(beta2 > param.Beta2Min))
        nb_racines_admissibles = nb_racines_admissibles + 1 ;
        racines_admissibles(nb_racines_admissibles) = racines(ii) ;
    end
end

% Recherche des points d'intersection de l'espace solution avec
% Beta2=Beta2min
[nbsols,racines] = RacinesReellesPolynomeDegre3(0,data.d2,data.d1,(data.d0t + data.dpi*PenteCom-param.Beta2Min)) ;
for ii=1:nbsols
    beta1 = data.c2*racines(ii)*racines(ii) + data.c1*racines(ii) + data.c0t + data.cpi*PenteCom;
    if ((abs(racines(ii)) < ThetaCible)&&(beta1 < param.Beta1Max)&&(beta1 > param.Beta1Min))
        nb_racines_admissibles = nb_racines_admissibles + 1 ;
        racines_admissibles(nb_racines_admissibles) = racines(ii) ;
    end
end

% Recherche des points d'intersection de l'espace solution avec
% Beta2=Beta2max
[nbsols,racines] = RacinesReellesPolynomeDegre3(0,data.d2,data.d1,(data.d0t + data.dpi*PenteCom-param.Beta2Max)) ;
for ii=1:nbsols
    beta1 = data.c2*racines(ii)*racines(ii) + data.c1*racines(ii) + data.c0t + data.cpi*PenteCom ;
    if ((abs(racines(ii)) < ThetaCible)&&(beta1 < param.Beta1Max)&&(beta1 > param.Beta1Min))
        nb_racines_admissibles = nb_racines_admissibles + 1 ;
        racines_admissibles(nb_racines_admissibles) = racines(ii) ;
    end
end

% --- definition des intervalles solution

% si pas de pt d'intersection : soit tous les points sont dedans, soit
% tous les points sont dehors ! On teste alors � Theta = 0
Beta1Test = data.cpi*PenteCom + data.c0t ;
Beta2Test = data.dpi*PenteCom + data.d0t ;
BetaTestOk = (Beta1Test > (param.Beta1Min+param.garde))&&(Beta1Test < (param.Beta1Max-param.garde))&&(Beta2Test > (param.Beta2Min+param.garde))&&(Beta2Test < (param.Beta2Max-param.garde)) ;

if ((nb_racines_admissibles>0)||(BetaTestOk))
    % alors la pente commandee est realisable et le sous-marin est
    % equilibrable
    valid = 1 ;
    PenteEq = PenteCom ;
else
    % La pente commandee n'est alors pas r�alisable et il faut trouver
    % comment modifier cette pente command�e

    % la pente est modifiee � la plus proche valeur admissible
    [valid,PenteEq] = RecherchePenteAdmissible(PenteCom,ThetaCible,param,data) ;

end



if (1==valid)
    data.c0 = data.c0t + data.cpi*PenteEq ;
    data.d0 = data.d0t + data.dpi*PenteEq ;
    % Enfin l'assiette d'equilibre est optimisee selon le critere
    [ThetaEq,Beta1Eq,Beta2Eq] = OptimisationAssietteEquilibre(PenteEq,ThetaCible,param,data) ;
else
    % Cas d�g�n�r�, sous-marin non equilibrable
    ThetaEq = 0 ;
    Beta1Eq = 0 ;
    Beta2Eq = 0 ;
end


