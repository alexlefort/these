function [c, intern] = tbx_sm_depth_controller(state,mission,intern,p,dt)
    
    % Interpolation des gains en fonction de la vitesse
    nVs    = length(p.gains_speeds);
    
    KZ     = tbx_math_interpolation(nVs, p.gains_speeds, p.gains_KZ     , state.U);
    IZ     = tbx_math_interpolation(nVs, p.gains_speeds, p.gains_IZ     , state.U);
    KPi    = tbx_math_interpolation(nVs, p.gains_speeds, p.gains_KPi    , state.U);
    KQ     = tbx_math_interpolation(nVs, p.gains_speeds, p.gains_KQ     , state.U);
    KTheta = tbx_math_interpolation(nVs, p.gains_speeds, p.gains_KTheta , state.U); 

    % Calcul Pente Commande
    Delta_Z     = (mission.ZCo - state.Z);
    PiCo        = KZ*Delta_Z;
    intern.Zint = IZ*Delta_Z*dt + intern.Zint + 0.05*intern.antiwp;
    PiCo        = PiCo + intern.Zint;
   
    PiCobrut      = PiCo;
    PiCo          = max(min(PiCo,p.sat_Pi_max),p.sat_Pi_min);
    intern.antiwp = PiCo - PiCobrut;

    % Calcul Pente
    Pi = state.Theta - atan(state.W/state.U);
    
    % Calcul de l ordre de barre brut
    c.Beta1 = (KPi+KTheta)*(PiCo - Pi) + 0*(PiCo - state.Theta) + KQ*(-state.Q);

    % Saturation et filtrage de l ordre de barre
    c.Beta1 = max(min(c.Beta1,p.sat_Beta),- p.sat_Beta);
    c.Beta1 = max(min(c.Beta1, intern.Beta + p.sat_Beta_d*dt), intern.Beta - p.sat_Beta_d*dt);
    intern.Beta = c.Beta1;
    c.Beta2 = 0.01;
    
    