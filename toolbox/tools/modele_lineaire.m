function [A,B,C,D] = modele_lineaire(Vs, param, nom, r)

axe.states = {'U';'V';'W';'P';'Q';'R';'Phi';'Theta';'Psi'};
axe.inputs = {'Alpha';'Beta1';'NT1'};
axe.mes    = {'U';'V';'W';'P';'Q';'R';'Phi';'Theta';'Psi'};

[Xeq, Ueq]         = calculer_equilibre(Vs, param, r);
[A_mdl,B_mdl,~,~]  = lineariser_modele(Xeq, Ueq, param, axe) ;

if (nom == 'GIR')
    %% Lateral : X = [V,P,R,Phi,Psi], U = ALpha
    A = A_mdl([2 4 6 7 9],[2 4 6 7 9]);
    B = B_mdl([2 4 6 7 9],1);
    C = eye(5);
    D = zeros(5,1);
end

if (nom == 'IMM')
    %% Longi : X = [W,Q,Theta], U = Beta1
    A = A_mdl([3 5 8],[3 5 8]);
    B = B_mdl([3 5 8],2);
    C = eye(3);
    D = zeros(3,1);

    ThetaPiQ = [0 0 1 ; -1/Vs 0 1 ; 0 1 0]; %% Transformee en X = [Theta, Pente, Q];

    A = ThetaPiQ*A*inv(ThetaPiQ);
    B = ThetaPiQ*B;
end
