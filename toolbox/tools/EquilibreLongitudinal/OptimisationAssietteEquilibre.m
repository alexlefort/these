function [ThetaEq,Beta1Eq,Beta2Eq] = OptimisationAssietteEquilibre(PenteEq,ThetaCible,param,data)
% Recherche de l'assiette optimale à PenteEq (supposée réalisable)

nb_racines_admissibles = 0 ;
racines_admissibles = [0,0,0,0,0,0,0,0] ;

nb_intervalles = 1 ; % au max 5, au min 1 car l'équilibre est réalisabl
bornes_inf_intervalles = [0,0,0,0,0] ;
bornes_sup_intervalles =  [0,0,0,0,0] ;

% Recherche des points d'intersection de l'espace solution avec
% Beta1=Beta1min
[nbsols,racines] = RacinesReellesPolynomeDegre3(0,data.c2,data.c1,(data.c0-param.Beta1Min)) ;
for ii=1:nbsols
    beta2 = data.d2*racines(ii)*racines(ii) + data.d1*racines(ii) + data.d0 ;
    if ((abs(racines(ii)) < ThetaCible)&&(beta2 < param.Beta2Max)&&(beta2 > param.Beta2Min))
        nb_racines_admissibles = nb_racines_admissibles + 1 ;
        racines_admissibles(nb_racines_admissibles) = racines(ii) ;
    end
end

% Recherche des points d'intersection de l'espace solution avec
% Beta1=Beta1max
[nbsols,racines] = RacinesReellesPolynomeDegre3(0,data.c2,data.c1,(data.c0-param.Beta1Max)) ;
for ii=1:nbsols
    beta2 = data.d2*racines(ii)*racines(ii) + data.d1*racines(ii) + data.d0 ;
    if ((abs(racines(ii)) < ThetaCible)&&(beta2 < param.Beta2Max)&&(beta2 > param.Beta2Min))
        nb_racines_admissibles = nb_racines_admissibles + 1 ;
        racines_admissibles(nb_racines_admissibles) = racines(ii) ;
    end
end

% Recherche des points d'intersection de l'espace solution avec
% Beta2=Beta2min
[nbsols,racines] = RacinesReellesPolynomeDegre3(0,data.d2,data.d1,(data.d0-param.Beta2Min)) ;
for ii=1:nbsols
    beta1 = data.c2*racines(ii)*racines(ii) + data.c1*racines(ii) + data.c0 ;
    if ((abs(racines(ii)) < ThetaCible)&&(beta1 < param.Beta1Max)&&(beta1 > param.Beta1Min))
        nb_racines_admissibles = nb_racines_admissibles + 1 ;
        racines_admissibles(nb_racines_admissibles) = racines(ii) ;
    end
end

% Recherche des points d'intersection de l'espace solution avec
% Beta2=Beta2max
[nbsols,racines] = RacinesReellesPolynomeDegre3(0,data.d2,data.d1,(data.d0-param.Beta2Max)) ;
for ii=1:nbsols
    beta1 = data.c2*racines(ii)*racines(ii) + data.c1*racines(ii) + data.c0 ;
    if ((abs(racines(ii)) < ThetaCible)&&(beta1 < param.Beta1Max)&&(beta1 > param.Beta1Min))
        nb_racines_admissibles = nb_racines_admissibles + 1 ;
        racines_admissibles(nb_racines_admissibles) = racines(ii) ;
    end
end


% --- definition des intervalles solution

if (0==nb_racines_admissibles)
    % necessairement, car l'equilibre est réalisable et aucun point
    % d'intersection avec les bordures du domaine n'a été trouvé
    nb_intervalles = 1 ;
    bornes_inf_intervalles(1) = -ThetaCible ;
    bornes_sup_intervalles(1) = ThetaCible ;
else

    pts_intersection = sort(racines_admissibles(1:nb_racines_admissibles)) ;

    % Test si le point -Thetamax appartient au domaine
    beta1 = data.c2*ThetaCible*ThetaCible - data.c1*ThetaCible + data.c0 ;
    beta2 = data.d2*ThetaCible*ThetaCible - data.d1*ThetaCible + data.d0 ;
    if ((beta1 < param.Beta1Max)&&(beta1 > param.Beta1Min)&&(beta2 < param.Beta2Max)&&(beta2 > param.Beta2Min))
        % si c'est le cas
        nb_intervalles = ceil((nb_racines_admissibles+1)/2) ;
        if (0==mod(nb_racines_admissibles,2))
            % si le nb de racines est pair : ThetaMax appartient au domaine
            bornes_inf_intervalles(1) = -ThetaCible ;
            bornes_sup_intervalles(nb_intervalles) = ThetaCible ;
            jj = 1 ;
            for ii=1:nb_racines_admissibles
                if (1==mod(ii,2)) % repartition des points d'intersection dans les bonnes cases
                    bornes_sup_intervalles(jj) = pts_intersection(ii) ;
                    jj= jj+1 ;
                else
                    bornes_inf_intervalles(jj) = pts_intersection(ii) ;
                end
            end

        else
            bornes_inf_intervalles(1) = -ThetaCible ;
            jj = 1 ;
            for ii=1:nb_racines_admissibles
                if (1==mod(ii,2)) % repartition des points d'intersection dans les bonnes cases
                    bornes_sup_intervalles(jj) = pts_intersection(ii) ;
                    jj= jj+1 ;
                else
                    bornes_inf_intervalles(jj) = pts_intersection(ii) ;
                end
            end

        end


    else
        % si -ThetaMax n'appartient pas au domaine
        nb_intervalles = floor((nb_racines_admissibles+1)/2) ;
        if (0==mod(nb_racines_admissibles,2))
            % si le nb de racines est pair : ThetaMax n'appartient pas au
            % domaine
            jj = 1 ;
            for ii=1:nb_racines_admissibles
                if (1==mod(ii,2)) % repartition des points d'intersection dans les bonnes cases
                    bornes_inf_intervalles(jj) = pts_intersection(ii) ;
                else
                    bornes_sup_intervalles(jj) = pts_intersection(ii) ;
                    jj= jj+1 ;
                end
            end

        else
            jj = 1 ;
            for ii=1:nb_racines_admissibles
                if (1==mod(ii,2)) % repartition des points d'intersection dans les bonnes cases
                    bornes_inf_intervalles(jj) = pts_intersection(ii) ;
                else
                    bornes_sup_intervalles(jj) = pts_intersection(ii) ;
                    jj= jj+1 ;
                end
            end
            bornes_sup_intervalles(nb_intervalles) = ThetaCible ;

        end

    end

end

% --- minimisation du critere J sur les intervalles solution

% Le critere a optimiser est */
% J = (b1/b10)^2 + (b2/b20)^2 + ((theta-pentecom)/eps0)^2+(theta/theta0)^2
% il peut se mettre sous la forme d'un polynome du quatrieme degre en fonction de l'assiette th
% J = e4.th^4+e3.th^3+e2.th^2+e1*th+e0
% les extremas sont donnees par la resolution de la derivee (dJ/dth) qui
% est un polynome du 3eme degre
e4 = (data.c2/param.b10)*(data.c2/param.b10)   + (data.d2/param.b20)*(data.d2/param.b20)   ;
e3 = (data.c2/param.b10)*(data.c1/param.b10)*2 + (data.d2/param.b20)*(data.d1/param.b20)*2 ;
e2 = (data.c1/param.b10)*(data.c1/param.b10)   + (data.d1/param.b20)*(data.d1/param.b20) + (data.c2/param.b10)*(data.c0/param.b10)*2 + (data.d2/param.b20)*(data.d0/param.b20)*2 + 1/(param.epsilon0*param.epsilon0) + 1/(param.theta0*param.theta0) ;
e1 = (data.c1/param.b10)*(data.c0/param.b10)*2 + (data.d1/param.b20)*(data.d0/param.b20)*2 - 2*PenteEq/(param.epsilon0*param.epsilon0) ;
[nbsols,racines_derivee] = RacinesReellesPolynomeDegre3(4*e4,3*e3,2*e2,e1) ;


J = 1e10 ; % l'infini par défaut
ThetaEq = 0 ;

for ii = 1:nb_intervalles
    % verification sur la borne inférieure de l'intervalle
   th_inf = bornes_inf_intervalles(ii) ;
   J_borne_inf = e4*th_inf^4 + e3*th_inf^3 + e2*th_inf^2 + e1*th_inf ;
   if (J_borne_inf < J)
       ThetaEq = th_inf ;
       J = J_borne_inf ;
   end
   % verification sur la borne supérieure de l'intervalle
   th_sup = bornes_sup_intervalles(ii) ;
   J_borne_sup = e4*th_sup^4 + e3*th_sup^3 + e2*th_sup^2 + e1*th_sup ;
   if (J_borne_sup < J)
       ThetaEq = th_sup ;
       J = J_borne_sup ;
   end


   % recherche des minimas dans l'intervalle
   for jj=1:nbsols
       thsol = racines_derivee(jj) ;
       derivee_seconde = 12*e4*thsol*thsol + 6*e3*thsol + 2*e2 ;
       if ((thsol > th_inf)&&(thsol < th_sup)&&(derivee_seconde >=0))
            J_sol = e4*thsol^4 + e3*thsol^3 + e2*thsol^2 + e1*thsol  ;
            if (J_sol < J)
                ThetaEq = thsol ;
                J = J_sol ;
            end
       end
   end
end

% Enfin, calcul des barres à l'equilibre
Beta1Eq = data.c2*ThetaEq*ThetaEq + data.c1*ThetaEq + data.c0 ;
Beta2Eq = data.d2*ThetaEq*ThetaEq + data.d1*ThetaEq + data.d0 ;
