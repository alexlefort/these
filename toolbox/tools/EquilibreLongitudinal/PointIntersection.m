function [valid,x0,y0] = PointIntersection(xm,ym,xu,yu,a,b,c)
% Calcule le point d'intersection des droites D1 et D2
% D1 : a.x + b.y + c = 0
% D2 : passant par (xm,ym) et de vect directeur (xu,yu) ;

delta = -a*xu - b*yu ;
if (0 == delta)
    % D1 // D2
    if (0==a*xm+b*ym+c)
        valid = 1 ;
        x0 = xm ;
        y0 = ym ;
    else
        valid = 0 ;
        x0 = 0 ;
        y0 = 0 ;
    end
else
    % il existe un unique point d'intersection
    valid = 1 ;
    x0 = (xu*c -b*(-xu*ym+yu*xm))/delta ;
    y0 = (yu*c + a*(-xu*ym+yu*xm))/delta ;
end
