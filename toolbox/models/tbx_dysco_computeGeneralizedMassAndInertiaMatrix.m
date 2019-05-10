function matriceGeneralisee = tbx_dysco_computeGeneralizedMassAndInertiaMatrix(C)
% TBX_DYSCO_COMPUTEGENERALIZEDMASSANDINERTIAMATRIX calcule la matrice
% g�n�ralis�e de masses et inerties ajout�es adimensionnalis�es par rho*Vol
%
%
% See also tbx_dysco_computeMassAndInertiaMatrix, tbx_dysco_computeAddedMassAndInertiaMatrix
%
% SIREHNA
% GJ
%==========================================================================
% SVN info
% SVN $Id$
% SVN $HeadURL$
%==========================================================================

matriceGeneralisee = ...
    tbx_dysco_computeMassAndInertiaMatrix(C) + ...
    tbx_dysco_computeAddedMassAndInertiaMatrix(C);
return;
