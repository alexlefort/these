function matriceInertie = tbx_dysco_computeMassAndInertiaMatrix(C)
% TBX_DYSCO_COMPUTEMASSANDINERTIAMATRIX evaluates the dimensionless
% mass and inertia matrix.
%
% \f[
% M_{\textrm{propre}} = \left[ {\begin{array}{*{20}c}
%    \mu  & 0 & 0 & 0 & {\mu Z_g } & { - \mu Y_g } \\
%    0 & \mu  & 0 & { - \mu Z_g } & 0 & {\mu X_g } \\
%    0 & 0 & \mu  & {\mu Y_g } & { - \mu X_g } & 0 \\
%    0 & { - \mu Z_g } & {\mu Y_g } & {L^2 \chi _1 } & 0 & { - L^2 \chi _{13} } \\
%    {\mu Z_g } & 0 & { - \mu X_g } & 0 & {L^2 \chi _2 } & 0 \\
%    { - \mu Y_g } & {\mu X_g } & 0 & { - L^2 \chi _{13} } & 0 & {L^2 \chi _3 } \\
% \end{array}} \right]
% \f]
%
% \note Cette matrice correspond à une matrice d'inertie dont tous les termes auraient
%       été divisés par \f[\rho.\textrm{Vol}\f] (Sorte d'adimensionnement pour se ramener au
%       coefficient de dépesée au lieu de la masse, avec \f[\textrm{Vol}\f] le volume de forme
%       du sous-marin
%
% See also tbx_dysco_computeAddedMassAndInertiaMatrix
%
% SIREHNA
% GJ
%==========================================================================
% SVN info
% SVN $Id$
% SVN $HeadURL$
%==========================================================================
L    = C.L;
Mu   = C.Mu;
X_1  = C.X_1;
X_2  = C.X_2;
X_3  = C.X_3;
X_13 = C.X_13;
Xg   = C.Xg;
Yg   = C.Yg;
Zg   = C.Zg;

matriceInertie = zeros(6,6);

matriceInertie(1,1) = +Mu;
matriceInertie(1,2) = 0;
matriceInertie(1,3) = 0;
matriceInertie(1,4) = 0;
matriceInertie(1,5) = +Mu * Zg;
matriceInertie(1,6) = -Mu * Yg;

matriceInertie(2,1) = 0;
matriceInertie(2,2) = +Mu;
matriceInertie(2,3) = 0;
matriceInertie(2,4) = -Mu * Zg;
matriceInertie(2,5) = 0;
matriceInertie(2,6) = +Mu * Xg;

matriceInertie(3,1) = 0;
matriceInertie(3,2) = 0;
matriceInertie(3,3) = +Mu;
matriceInertie(3,4) = +Mu * Yg;
matriceInertie(3,5) = -Mu * Xg;
matriceInertie(3,6) = 0;

matriceInertie(4,1) = 0;
matriceInertie(4,2) = -Mu  * Zg;
matriceInertie(4,3) = +Mu  * Yg;
matriceInertie(4,4) = +L*L * X_1;
matriceInertie(4,5) = 0;
matriceInertie(4,6) = -L*L * X_13;

matriceInertie(5,1) = +Mu  * Zg;
matriceInertie(5,2) = 0;
matriceInertie(5,3) = -Mu  * Xg;
matriceInertie(5,4) = 0;
matriceInertie(5,5) = +L*L * X_2;
matriceInertie(5,6) = 0;

matriceInertie(6,1) = -Mu  * Yg;
matriceInertie(6,2) = +Mu  * Xg;
matriceInertie(6,3) = 0;
matriceInertie(6,4) = -L*L * X_13;
matriceInertie(6,5) = 0;
matriceInertie(6,6) = L*L  * X_3;
return;
