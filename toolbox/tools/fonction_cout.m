function cout = fonction_cout(xu)

% fonction de cout pour la recherche d'equilibre
% xu = [U, V, W, P, Q, R, Phi, Theta, Alpha, Beta1, Beta2, NT1]'

global p;
global xobj;

[x,u] = init_etats;

xf = fields(x);
uf = fields(u);

for jj=1:8
    x.(xf{jj}) = xu(jj);
end

for jj=1:4
    u.(uf{jj}) = xu(8 + jj);
end

[,xpt] = equations_modele(x,u,p,0.1) ;

x1 = [x.U; x.R; xpt.Z; x.Theta] ; % U, R, Zpt, Theta
matadim1 = diag([1/xobj(1) p.L/xobj(1) 2/xobj(1) 1]) ;

x2pt = [xpt.U; xpt.V; xpt.W; xpt.P; xpt.Q; xpt.R; xpt.Phi; xpt.Theta] ;

matadim2 = diag([1/xobj(1) 1/xobj(1) 1/xobj(1) p.L/xobj(1) p.L/xobj(1) p.L/xobj(1) 1 1]) ;

k = 1*p.L/xobj(1) ;
cout = (matadim1*(x1-xobj))' * (matadim1*(x1-xobj)) + k^2 * (matadim2*x2pt)' * (matadim2*x2pt) ;
