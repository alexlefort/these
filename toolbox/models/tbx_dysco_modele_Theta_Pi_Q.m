function [A,B] = tbx_dysco_modele_Theta_Pi_Q(param,Vs,g)
% tbx_dysco_modele_Theta_W_Q returns the state space matrix of a linear model [Theta,Pi,Q,beta1,beta2]
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

[A1,B1] = tbx_dysco_modele_Theta_W_Q(param,Vs,g);
% transform:
%   X = [Theta ; Pi ; Q], u= [beta1 ; beta2]  Pi = theta - W/Vs
Mpass = [1 0 0 ; 1 -1/Vs 0 ; 0 0 1] ;
A = Mpass*A1*inv(Mpass) ;
B = Mpass*B1 ;