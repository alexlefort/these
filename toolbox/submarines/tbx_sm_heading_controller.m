function [c, intern] = tbx_sm_heading_controller(state,mission,intern,p,dt)

    % Interpolation des gains en fonction de la vitesse
    nVs  = length(p.gains_speeds);
    KPSI = tbx_math_interpolation(nVs,p.gains_speeds,p.gains_KPSI,state.U);
    KR   = tbx_math_interpolation(nVs,p.gains_speeds,p.gains_KR  ,state.U);
    IR   = tbx_math_interpolation(nVs,p.gains_speeds,p.gains_IR  ,state.U);
    
    % Calcul Taux de Giration Commande
    RCo  = KPSI*(mod(mission.PsiCo - state.Psi + pi , 2*pi) - pi);
    RCo  = max(min(RCo,p.sat_R_max),-p.sat_R_max);
    
    % Calcul de l ordre de barre brut
    Delta_R = RCo - state.R;
    c.Alpha = KR*Delta_R;
    intern.Rint = IR*Delta_R*dt + intern.Rint + intern.antiwp;
    c.Alpha = c.Alpha + intern.Rint;
    
    % Saturation et filtrage de l ordre de barre
    alpha_brut = c.Alpha;
    c.Alpha = max(min(c.Alpha,p.sat_Alpha),- p.sat_Alpha);
    c.Alpha = max(min(c.Alpha, intern.Alpha + p.sat_Alpha_d*dt), intern.Alpha - p.sat_Alpha_d*dt);
    intern.antiwp = c.Alpha - alpha_brut;
    intern.Alpha = c.Alpha;
    