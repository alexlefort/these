function [intersection,delta_pente] = IncrementPente(beta1,beta2,param,data)
% calcule la variation de pente d'equilibre pour ramener le point
% (beta1,beta2) dans le domaine admissible

intersection = 0 ;
delta_pente = 1e10 ; % initialisation à l' "infini"

lambda = data.cpi*data.cpi + data.dpi*data.dpi ;
if (lambda > 0) % si le pb est sensible à la pente d'équilibre
    % calcule point d'intersection avec beta1 = beta1min+garde
    a = 1 ;
    b = 0 ;
    c = -(param.Beta1Min+param.garde) ;
    [valid,x0,y0] = PointIntersection(beta1,beta2,data.cpi,data.dpi,a,b,c) ;
    if valid && (y0 < (param.Beta2Max-param.garde))&&(y0 > (param.Beta2Min+param.garde))
        % le point d'intesection est sur le segment de Beta2
        intersection = 1 ;
        dp = ((x0-beta1)*data.cpi + (y0-beta2)*data.dpi)/lambda ;
        if (abs(dp)<abs(delta_pente)) % cette variation est plus interessante
            delta_pente = dp ;
        end
    end

    % calcule point d'intersection avec beta1 = beta1max-garde
    a = 1 ;
    b = 0 ;
    c = -(param.Beta1Max-param.garde) ;
    [valid,x0,y0] = PointIntersection(beta1,beta2,data.cpi,data.dpi,a,b,c) ;
    if valid && (y0 < (param.Beta2Max-param.garde))&&(y0 > (param.Beta2Min+param.garde))
        % le point d'intesection est sur le segment de Beta2
        intersection = 1 ;
        dp = ((x0-beta1)*data.cpi + (y0-beta2)*data.dpi)/lambda ;
        if (abs(dp)<abs(delta_pente)) % cette variation est plus interessante
            delta_pente = dp ;
        end
    end

    % calcule point d'intersection avec beta2 = beta2min+garde
    a = 0 ;
    b = 1 ;
    c = -(param.Beta2Min+param.garde) ;
    [valid,x0,y0] = PointIntersection(beta1,beta2,data.cpi,data.dpi,a,b,c) ;
    if valid && (x0 < (param.Beta1Max-param.garde))&&(x0 > (param.Beta1Min+param.garde))
        % le point d'intesection est sur le segment de Beta1
        intersection = 1 ;
        dp = ((x0-beta1)*data.cpi + (y0-beta2)*data.dpi)/lambda ;
        if (abs(dp)<abs(delta_pente)) % cette variation est plus interessante
            delta_pente = dp ;
        end
    end

    % calcule point d'intersection avec beta2 = beta2max-garde
    a = 0 ;
    b = 1 ;
    c = -(param.Beta2Max-param.garde) ;
    [valid,x0,y0] = PointIntersection(beta1,beta2,data.cpi,data.dpi,a,b,c) ;
    if valid && (x0 < (param.Beta1Max-param.garde))&&(x0 > (param.Beta1Min+param.garde))
        % le point d'intesection est sur le segment de Beta1
        intersection = 1 ;
        dp = ((x0-beta1)*data.cpi + (y0-beta2)*data.dpi)/lambda ;
        if (abs(dp)<abs(delta_pente)) % cette variation est plus interessante
            delta_pente = dp ;
        end
    end

end