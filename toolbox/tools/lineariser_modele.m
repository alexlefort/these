function [A,B,C,D] = lineariser_modele(Xeq, Ueq, param, axe) ;

deg = pi/180;
dd  = 1e-6;
dd2 = 1e6;

ns = numel(axe.states);
ni = numel(axe.inputs);
nm = numel(axe.mes   );

% Calcul des gradients
A = zeros(ns, ns) ;
B = zeros(ns, ni) ;
C = zeros(nm, ns) ;
D = zeros(nm, ni) ;

X0 = Xeq;
U0 = Ueq;
X  = X0;
U  = U0;

% Calcul au point d'equilibre
[y0,DX0] = equations_modele(Xeq, Ueq, param, 0.1);

for ii=1:ns
    state          = axe.states{ii};
    X.(state)      = X.(state) + dd;
    [y,DX.(state)] = equations_modele(X, U, param, 0.1);
    X              = X0;
end

for ii=1:ni
    input          = axe.inputs{ii};
    U.(input)     = U.(input) + dd;
    [y,DX.(input)] = equations_modele(X, U, param, 0.1);
    U             = U0;
end

%% Matrice d'etat

for ii=1:ns
    for jj=1:ns
        state_ii = axe.states{ii};
        state_jj = axe.states{jj};
        A(ii,jj) = (DX.(state_jj).(state_ii) - DX0.(state_ii))*dd2;
    end
end

%% Matrice de commande

for ii=1:ns
    for jj=1:ni
        state = axe.states{ii,1};
        input = axe.inputs{jj,1};
        B(ii,jj) = (DX.(input).(state) - DX0.(state))*dd2;
    end
end

%% Matrice d'observation

for ii=1:nm
    for jj=1:ns
    mes   = axe.mes(ii);
    state = axe.states(jj);
        C(ii,jj) = strcmp(mes,state);
    end
end
