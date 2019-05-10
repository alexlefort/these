function [A,B] = tbx_dysco_modele_Phi_Psi_V_P_R(param,Vs)
% tbx_dysco_modele_Phi_Psi_V_P_R returns the state space matrix of a linear model [Phi,Psi,V,P,R,alpha]
%
% Inputs:
%   - param: Structure containing all coefficients from man model.
%            typically as obtained with tbx_dysco_readNavData
%   - Vs: velocity at which the model is linearized.
%
% Outputs:
%   - A,B : matrices X dot = A*X + B*U
%           where X=[ Phi; Psi; V ; P ; R] and U=Alpha
%
% Requirements:
% - The model must be in cross configuration (2 actuators alpha, beta1)
% - The case of parameters must be normalized (See tbx_dysco_normalizeNavDataCase)
%
% Note that this is a simplified model.
% Linearisation point is U=Vs W=V=0 P=Q=R=0 Phi=Theta=0
% Higher order terms in manoeuvrability model are neglected
%
%==========================================================================
% SVN info
% SVN $Id$
% SVN $HeadURL$
%==========================================================================

% matrices de masse et masses ajoutees
matriceGeneralisee = tbx_dysco_computeGeneralizedMassAndInertiaMatrix(param);

% Projection sur roulis-lacet  V-P-R
Mlat2 = matriceGeneralisee([2 4 6],[2 4 6]);
Minvlat2 = [eye(2) zeros(2,3) ; zeros(3,2) inv(Mlat2)] ;

% creation d'un mdl lineaire roulis et lacet : X dot = A*X + B*U
%   X = [Phi Psi V P R], U=[ alpha ]
Ainert2 = [  0 0 0 1 0 ; 0 0 0 0 1 ; 0 0 0 0 -param.Mu*Vs ; 0 0 0 0 param.Mu*param.Zg*Vs ; 0 0 0 0 -param.Mu*param.Xg*Vs  ] ;
Ahydrostat2 = [ zeros(2,5) ; (param.Mu-1)*param.g 0 0 0 0 ; -param.Mu*param.Zg*param.g 0 0 0 0 ; param.Mu*param.Xg*param.g 0 0 0 0] ;
Aetat2 = Vs/(2*param.CRemp)*[ zeros(2,5) ; 0 0 param.CyV/param.L param.CyP param.CyR  ;  0 0 param.ClV param.ClP*param.L param.ClR*param.L ; 0 0 param.CnV param.CnP*param.L param.CnR*param.L ] ;
A = Minvlat2*(Ainert2+Aetat2+Ahydrostat2) ;
B = Minvlat2*Vs*Vs/(2*param.CRemp)*[0 ; 0 ; param.CyAL/param.L ; param.ClAL ; param.CnAL] ;
