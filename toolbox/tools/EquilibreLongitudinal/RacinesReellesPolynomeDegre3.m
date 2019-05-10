function [nbsols,racines] = RacinesReellesPolynomeDegre3(a,b,c,d)
% Resolution dans IR du polynome a.x^3 + b.x^2 + c.x + d = 0
% nbsols : indique le nombre de solutions valides retournees par la fonction
% (ex : si nbsols==2, l'ensemble des solutions est {racines(1), racines(2)}
%       si nbsols==0, l'ensemble des solutions est vide)
% Les solutions valides sont donnees dans l'ordre croissant
% On suppose (a,b,c) different de (0,0,0)

racines = [0,0,0] ;

if (0==a)
    if (0==b)
        if (0==c)
            % le polynome est de degre 0 : cas non traité
            nbsols = 0 ;
        else
            % le polynome est de degre 1
            nbsols = 1 ;
            racines(1) = -d/c ;
        end
    else
        % le polynome est de degré 2
        delta = c*c - 4*b*d ;
        if (delta < 0)
            nbsols = 0 ;
        elseif (0==delta)
            nbsols = 1 ;
            racines(1) = -c/(2*b) ;
        else
            nbsols = 2 ;
            racines(1) = (-c-sign(b)*sqrt(delta))/(2*b) ;
            racines(2) = (-c+sign(b)*sqrt(delta))/(2*b) ;
        end
    end
else
    % le polynome est de degré 3
    p = -b*b/(3*a*a) + c/a ;
    q = b/(27*a)*(2*b*b/(a*a)-9*c/a)+d/a ;
    delta = q*q+4*p*p*p/27 ;
    if (delta > 0)
        nbsols = 1 ;
        sign1 = sign(-q+sqrt(delta)) ;
        sign2 = sign(-q-sqrt(delta)) ;
        racines(1) = -b/(3*a) + sign1*power(abs(-q/2+sqrt(delta)/2),1/3)+sign2*power(abs(-q/2-sqrt(delta)/2),1/3) ;
    elseif (delta < 0)
        % necessairement p est different de 0
        nbsols = 3 ;
        racines(1) = -b/(3*a) + 2*sqrt(-p/3)*cos(acos(-q/2*sqrt(-27/(p*p*p)))/3);
        racines(2) = -b/(3*a) + 2*sqrt(-p/3)*cos(acos(-q/2*sqrt(-27/(p*p*p)))/3 + 2*pi/3);
        racines(3) = -b/(3*a) + 2*sqrt(-p/3)*cos(acos(-q/2*sqrt(-27/(p*p*p)))/3 + 4*pi/3);
        racines = sort(racines) ;
    else
        if (0==p)
            nbsols = 1 ;
            racines(1) = -b/(3*a) ;
        else
            nbsols = 2 ;
            if (q*p < 0)
                racines(1) = -b/(3*a) + 3*q/p;
                racines(2) = -b/(3*a) - 3*q/(2*p);
            else
                racines(1) = -b/(3*a) - 3*q/(2*p);
                racines(2) = -b/(3*a) + 3*q/p;
            end
        end
    end
end
