function [y,xpt] = equations_modele(x,u,p,dt)

% Calcul du modele DYSCO pour
% x = [U, V, W, P, Q, R, X, Y, Z, Phi, Theta, Psi, fxp, fyp, fzp, mxp, myp, mzp]'
% u = [Alpha, Beta1, Beta2, NT1]'
% p : structure avec les parametres fichiers .nav

U      = x.U     ;
V      = x.V     ;
W      = x.W     ;
P      = x.P     ;
Q      = x.Q     ;
R      = x.R     ;
Phi    = x.Phi   ;
Theta  = x.Theta ;
Psi    = x.Psi   ;
fxp    = x.fxp   ;
fzp    = x.fzp   ;
myp    = x.myp   ;

Alpha = u.Alpha  ;
Beta1 = u.Beta1  ;
Beta2 = 0        ;
NT1   = u.NT1    ;

% Matrice masses + masses ajoutees
M = [          (p.Mu+p.Mu_1)                        0                  p.Mu_13                           0  (p.Mu*p.Zg-p.L*p.Nu_15)                  -p.Mu*p.Yg ;...
                           0            (p.Mu+p.Mu_2)                        0    (-p.Mu*p.Zg-p.L*p.Nu_24)                        0     (p.Mu*p.Xg-p.L*p.Nu_26) ;...
                     p.Mu_13                        0            (p.Mu+p.Mu_3)                   p.Mu*p.Yg (-p.Mu*p.Xg-p.L*p.Nu_35)                           0 ;...
                           0 (-p.Mu*p.Zg-p.L*p.Nu_24)                p.Mu*p.Yg    p.L^2*(p.X_1+p.Lambda_1)                        0 -p.L^2*(p.X_13+p.Lambda_13) ;...
     (p.Mu*p.Zg-p.L*p.Nu_15)                        0 (-p.Mu*p.Xg-p.L*p.Nu_35)                           0 p.L^2*(p.X_2+p.Lambda_2)                           0 ;...
                  -p.Mu*p.Yg  (p.Mu*p.Xg-p.L*p.Nu_26)                        0 -p.L^2*(p.X_13+p.Lambda_13)                        0    p.L^2*(p.X_3+p.Lambda_3) ] ;

M = [ M zeros(6) ; zeros(6) eye(6)] ; % ajout des 6 autres dimensions

% Effort Hydrostatique
FHydroStatique = [ (1-p.Mu)*p.g*sin(Theta) ;...
                   (p.Mu-1)*p.g*cos(Theta)*sin(Phi) ;...
                   (p.Mu-1)*p.g*cos(Theta)*cos(Phi) ;...
                   -p.Mu*p.Zg*p.g*cos(Theta)*sin(Phi) + p.Mu*p.Yg*p.g*cos(Theta)*cos(Phi) ;...
                   -p.Mu*p.Xg*p.g*cos(Theta)*cos(Phi) - p.Mu*p.Zg*p.g*sin(Theta) ;...
                   p.Mu*p.Xg*p.g*cos(Theta)*sin(Phi) + p.Mu*p.Yg*p.g*sin(Theta)
                   zeros(6,1)] ;

% Effort couplage Inertiel
FCouplageInertiel = [ p.Mu*(V*R-W*Q + p.Xg*(Q^2+R^2)-p.Yg*P*Q-p.Zg*P*R) ;...
                      p.Mu*(W*P-U*R-p.Xg*P*Q+p.Yg*(P^2+R^2)-p.Zg*Q*R) ;...
                      p.Mu*(U*Q-V*P-p.Xg*P*R-p.Yg*Q*R+p.Zg*(P^2+Q^2)) ;...
                      p.L^2*((p.X_2-p.X_3)*Q*R + p.X_13*P*Q) + p.Mu*(p.Yg*(U*Q-V*P)+p.Zg*(U*R-W*P)) ;...
                      p.L^2*((p.X_3-p.X_1)*P*R + (p.X_13+p.Lambda_13)*(R^2-P^2)) + p.Mu*(-p.Xg*(U*Q-V*P)+p.Zg*(V*R-W*Q)) ;...
                      p.L^2*((p.X_1-p.X_2)*P*Q - p.X_13*Q*R) + p.Mu*(p.Xg*(W*P-U*R)+p.Yg*(W*Q-V*R)) ;...
                      U*cos(Psi)*cos(Theta)+V*(cos(Psi)*sin(Theta)*sin(Phi)-sin(Psi)*cos(Phi))+W*(sin(Psi)*sin(Phi)+cos(Psi)*sin(Theta)*cos(Phi)) ;...
                      U*sin(Psi)*cos(Theta) + V*(cos(Psi)*cos(Phi)+sin(Psi)*sin(Theta)*sin(Phi)) + W*(sin(Psi)*sin(Theta)*cos(Phi)-cos(Psi)*sin(Phi)) ;...
                      -U*sin(Theta)+V*cos(Theta)*sin(Phi) + W*cos(Theta)*cos(Phi) ;...
                      P+(Q*sin(Phi)+R*cos(Phi))*tan(Theta) ;...
                      Q*cos(Phi)-R*sin(Phi) ;...
                      (R*cos(Phi)+Q*sin(Phi))/cos(Theta) ] ;

% Effort hydrodynamique
Vs = sqrt(U^2+V^2+W^2) ;

ured = U/Vs ;
vred = V/Vs ;
wred = W/Vs ;
pred = p.L*P/Vs;
qred = p.L*Q/Vs ;
rred = p.L*R/Vs ;

Cx = p.Cx0 ;
Cy = p.Cy0 + p.CyV*vred + p.CyR*rred + p.CyP*pred + p.CyAL*Alpha ;
Cz = p.Cz0 + p.CzW*wred + p.CzQ*qred + p.CzB1*Beta1 + p.CzB2*Beta2 ;
Cl = p.Cl0 + p.ClV*vred + p.ClR*rred + p.ClP*pred + p.ClAL*Alpha ;
Cm = p.Cm0 + p.CmW*wred + p.CmQ*qred + p.CmB1*Beta1 + p.CmB2*Beta2 ;
Cn = p.Cn0 + p.CnV*vred + p.CnR*rred + p.CnAL*Alpha;
FHydrodynamique = Vs^2/(2*p.CRemp)*[ Cx/p.L ;...
                                     Cy/p.L ;...
                                     Cz/p.L ;...
                                     Cl  ; ...
                                     Cm  ; ...
                                     Cn  ; 0 ; 0 ; 0 ; 0 ; 0 ; 0];


FPerturbateur = [ fxp ; 0 ; fzp ; 0 ; myp ; 0 ; zeros(6,1)];

% Effort de propulsion

Js = Vs/(NT1*p.D) ;
Js_etoile = Js*(1+p.Aj*sqrt((vred+p.XAR*rred)^2+(wred-p.XAR*qred)^2) + p.Bj*sqrt((vred+p.XAR*rred)^2+(wred-p.XAR*qred)^2)) ;
Vs_etoile = Js_etoile*NT1*p.D ;

if (abs(Js_etoile) <= 1)
    Teffred = NT1^2*p.D^4*interpolation(10,p.J,p.KV,Js_etoile)/(p.CRemp*p.A*p.L) ;
    Qeffred = -NT1^2*p.D^5*interpolation(10,p.J,p.KQ,Js_etoile)/(p.CRemp*p.A*p.L) ;
else
    Teffred = p.D^2*Vs_etoile^2*interpolation(10,p.unSurJ,p.CV,1/Js_etoile)/(p.CRemp*p.A*p.L)  ;
    Qeffred = p.D^3*Vs_etoile^2*interpolation(10,p.unSurJ,p.CQ,1/Js_etoile)/(p.CRemp*p.A*p.L)  ;
end

FPropulsion = [Teffred 0 0 Qeffred 0 0 0 0 0 0 0 0 ]' ;

% Equation d'evolution
xpt_vect = inv(M)*(FHydroStatique + FCouplageInertiel + FHydrodynamique + FPropulsion + FPerturbateur) ;

xpt.U      = xpt_vect(1)  ;
xpt.V      = xpt_vect(2)  ;
xpt.W      = xpt_vect(3)  ;
xpt.P      = xpt_vect(4)  ;
xpt.Q      = xpt_vect(5)  ;
xpt.R      = xpt_vect(6)  ;
xpt.Phi    = xpt_vect(10) ;
xpt.Theta  = xpt_vect(11) ;
xpt.Psi    = xpt_vect(12) ;
xpt.Z      = xpt_vect(9)  ;
xpt.fxp    = 0.0          ;
xpt.fzp    = 0.0          ;
xpt.myp    = 0.0          ;

y.U        = x.U     + dt * xpt.U     ;
y.V        = x.V     + dt * xpt.V     ;
y.W        = x.W     + dt * xpt.W     ;
y.P        = x.P     + dt * xpt.P     ;
y.Q        = x.Q     + dt * xpt.Q     ;
y.R        = x.R     + dt * xpt.R     ;
y.Phi      = x.Phi   + dt * xpt.Phi   ;
y.Theta    = x.Theta + dt * xpt.Theta ;
y.Psi      = x.Psi   + dt * xpt.Psi   ;
y.Z        = x.Z     + dt * xpt.Z     ;
y.fxp      = x.fxp   + dt * xpt.fxp   ;
y.fzp      = x.fzp   + dt * xpt.fzp   ;
y.myp      = x.myp   + dt * xpt.myp   ;
