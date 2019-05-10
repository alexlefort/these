function matriceAjoutee = tbx_dysco_computeAddedMassAndInertiaMatrix(C)
% TBX_DYSCO_COMPUTEADDEDMASSANDINERTIAMATRIX calcule la matrice de masses et
% inerties ajoutées adimensionnalisées par rho*Vol
%
% \f[
% \left[ {M_{\textrm{ajout\'ee}} } \right] = \left[ {\begin{array}{*{20}c}
%    {\mu _1 } & 0 & {\mu _{13} } & 0 & { - L\nu _{15} } & 0 \\
%    0 & {\mu _2 } & 0 & { - L\nu _{24} } & 0 & { - L\nu _{26} } \\
%    {\mu _{13} } & 0 & {\mu _3 } & 0 & { - L\nu _{35} } & 0 \\
%    0 & { - L\nu _{24} } & 0 & {L^2 \lambda _1 } & 0 & { - L^2 \lambda _{13} } \\
%    { - L\nu _{15} } & 0 & { - L\nu _{35} } & 0 & {L^2 \lambda _2 } & 0 \\
%    0 & { - L\nu _{26} } & 0 & { - L^2 \lambda _{13} } & 0 & {L^2 \lambda _3 } \\
% \end{array}} \right]
% \f]
%
% \param[out] matriceAjoutee Matrice de masse
% \param[in] C Paramètres d'inertie
%
% See also tbx_dysco_computeMassAndInertiaMatrix
%
% SIREHNA
% GJ
%==========================================================================
% SVN info
% SVN $Id$
% SVN $HeadURL$
%==========================================================================
L         = C.L;
Mu_1      = C.Mu_1;
Mu_13     = C.Mu_13;
Mu_2      = C.Mu_2;
Mu_3      = C.Mu_3;
Nu_15     = C.Nu_15;
Nu_24     = C.Nu_24;
Nu_26     = C.Nu_26;
Nu_35     = C.Nu_35;
Lambda_1  = C.Lambda_1;
Lambda_13 = C.Lambda_13;
Lambda_2  = C.Lambda_2;
Lambda_3  = C.Lambda_3;

matriceAjoutee = zeros(6,6);

matriceAjoutee(1,1) = Mu_1;
matriceAjoutee(1,2) = 0;
matriceAjoutee(1,3) = Mu_13;
matriceAjoutee(1,4) = 0;
matriceAjoutee(1,5) = - L * Nu_15;
matriceAjoutee(1,6) = 0;

matriceAjoutee(2,1) = 0;
matriceAjoutee(2,2) = Mu_2;
matriceAjoutee(2,3) = 0;
matriceAjoutee(2,4) = - L * Nu_24;
matriceAjoutee(2,5) = 0;
matriceAjoutee(2,6) = - L * Nu_26;

matriceAjoutee(3,1) = Mu_13;
matriceAjoutee(3,2) = 0;
matriceAjoutee(3,3) = Mu_3;
matriceAjoutee(3,4) = 0;
matriceAjoutee(3,5) = - L * Nu_35;
matriceAjoutee(3,6) = 0;

matriceAjoutee(4,1) = 0;
matriceAjoutee(4,2) = - L * Nu_24;
matriceAjoutee(4,3) = 0;
matriceAjoutee(4,4) = L * L*Lambda_1;
matriceAjoutee(4,5) = 0;
matriceAjoutee(4,6) = -L*L*Lambda_13;

matriceAjoutee(5,1) = - L * Nu_15;
matriceAjoutee(5,2) = 0;
matriceAjoutee(5,3) = - L * Nu_35;
matriceAjoutee(5,4) = 0;
matriceAjoutee(5,5) = L*L * Lambda_2;
matriceAjoutee(5,6) = 0;

matriceAjoutee(6,1) = 0;
matriceAjoutee(6,2) = - L * Nu_26;
matriceAjoutee(6,3) = 0;
matriceAjoutee(6,4) = -L*L*Lambda_13;
matriceAjoutee(6,5) = 0;
matriceAjoutee(6,6) = L*L*Lambda_3;
return;
