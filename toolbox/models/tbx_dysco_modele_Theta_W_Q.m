function [A,B] = tbx_dysco_modele_Theta_W_Q(param,Vs,g)
% tbx_dysco_modele_Theta_W_Q returns the state space matrix of a linear model [Theta,W,Q,beta1,beta2]
%
% Inputs:
%   - param: Structure containing all coefficients from man model.
%            typically as obtained with tbx_dysco_readNavData
%   - Vs: velocity at which the model is linearized.
%   - g (optional): acceleration of gravity (default to 9.81 m/s2)
%
% Outputs:
%   - A,B : matrices X dot = A*X + B*U
%           where X=[ Theta ; W ; Q ] and U=[ beta1 ; beta2 ]
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

if nargin<3
    if isfield(param,'g')
        g=param.g;
    else
        g=9.81;
    end
end

% matrices de masse et masses ajoutées
matriceGeneralisee = tbx_dysco_computeGeneralizedMassAndInertiaMatrix(param);

% projection sur les axes W et Q    
Mlong = matriceGeneralisee([3 5],[3 5]);           
Minvlong = inv([1 0 0 ; zeros(2,1) Mlong]) ;

% modele X dot = A*X+B*U avec
%   X = [Theta ; W ; Q], U= [beta1 ; beta2]
Ainert = [  0 0 1 ; 0 0 param.Mu*Vs ; 0 0 -param.Mu*param.Xg*Vs] ;
Ahydrostat = [ 0 0 0 ; 0 0 0 ; -param.Mu*param.Zg*g 0 0] ;
Aetat = Vs/(2*param.CRemp)*[ 0 0 0 ; 0 param.CzW/param.L param.CzQ  ;  0 param.CmW  param.CmQ*param.L] ;
A = Minvlong*(Ainert+Aetat+Ahydrostat) ;
B = Minvlong*Vs*Vs/(2*param.CRemp)*[0 0 ; param.CzB1/param.L param.CzB2/param.L; param.CmB1 param.CmB2] ;

