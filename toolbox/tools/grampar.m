%% Calcul du grammien de commandabilité partiel :
%% Temps de reponse d'un systeme regle par LQ

%% Entrees :
%% sys : systeme a asservir
%% Tc  : temps de reponse caracteristique
%% dt  : pas de temps pour le calcul

function result = grampar(sys,Tc,dt)

A = sys.a;
B = sys.b;
C = sys.c;

vect_t = dt:dt:Tc;
size_t = size(vect_t,2);

result = expm(A'*0)*C'*C*expm(A*0)*dt;

for ii=1:size_t
    t = vect_t(ii);
    result = result + expm(A'*t)*C'*C*expm(A*t)*dt;
end
