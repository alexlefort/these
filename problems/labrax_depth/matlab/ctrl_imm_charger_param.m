%% p de la loi de commande Immersion

function p = ctrl_imm_charger_param

deg = pi/180;

%% Saturations de barres

p.sat_Beta    =  20*deg ;
p.sat_Beta_d  =  12*deg ;

%% Saturations de commande

p.sat_Z_min   =    0     ;
p.sat_Z_max   =  100     ;
p.sat_Pi_min  =  -15*deg ;
p.sat_Pi_max  =   15*deg ;

%% Anti-windup

p.gain_antiwp = 1.0;
%% Parametres du module d equilibre

%% Reglage de la methode de reglage

% modelisation d'un retard pur pour le reglage des gains
% retard de 2s/sqrt(echelle) modelise a l ordre 4

[p.num_pade,p.den_pade] = pade(0.5,4) ;

%% Points de réglages

%% Vs_Vect    : tabulation en vitesse
%% KZmin_Vect : gain KZ max considere : on augmente le gain jusqu a ce le depassement disparaisse.
%% Trep_Vect  : temps de reponse moyen voulu pour les etats de la boucle interne.
p.Tvect = 0:0.01:100;
	
%% Points de reglages

p.Vs_Vect     = [   0.3   0.5  1.0   1.5   2.0   2.5   3.0   4.0   5.0];
p.Trep_Vect   = [  10.0   5.0  3.0   2.7   2.5   2.3   2.1   1.9   1.7];

p.nVs         = length(p.Vs_Vect);
