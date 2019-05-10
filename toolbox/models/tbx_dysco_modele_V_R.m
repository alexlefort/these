function [A,B] = tbx_dysco_modele_V_R(param,Vs)
% tbx_dysco_modele_V_R returns the state space matrix of a linear model [V,R,alpha]
%
% Inputs:
%   - param: Structure containing all coefficients from man model.
%            typically as obtained with tbx_dysco_readNavData
%   - Vs: velocity at which the model is linearized.
%
% Outputs:
%   - A,B : matrices X dot = A*X + B*U
%           where X=[ V ; R ] and U=Alpha
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

% projection sur les axes V et R    
Mlat = matriceGeneralisee([2 6],[2 6]);           
Minvlat = inv(Mlat) ;

% modele X dot = A*X+B*U avec
%   X = [V ; R], U= [ alpha ]
Ainert = [  0 -param.Mu*Vs ; 0 -param.Mu*param.Xg*Vs] ;
Aetat = Vs/(2*param.CRemp)*[ param.CyV/param.L param.CyR  ;  param.CnV  param.CnR*param.L] ;
A = Minvlat*(Ainert+Aetat) ;
B = Minvlat*Vs*Vs/(2*param.CRemp)*[ param.CyAL/param.L ; param.CnAL] ;

