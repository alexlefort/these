%% A disposition :

%% Amat
%% Bmat
%% seq
%% ueq
%% s
%% u

param = ChargerModeleSM();


%% Comparaison des modèles formels et numériques

deg = pi/180;

Vsvect  = [ 0.5 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0] ; %% Table des vitesses
nVs = length(Vsvect) ;

ALL.states = {'U';'V';'W';'P';'Q';'R';'X';'Y';'Z';'Phi';'Theta';'Psi'}  ;  ALL.nb_states = size(ALL.states,1);
ALL.inputs = {'Alpha';'Beta1';'Beta2';'NT1'}                            ;  ALL.nb_inputs = size(ALL.inputs,1);
	
GIR.states = {'V';'P';'R';'Phi';'Psi'} ;  GIR.nb_states = size(GIR.states,1);
GIR.inputs = {'Alpha'}                 ;  GIR.nb_inputs = size(GIR.inputs,1);

Amat_gir = Amat([2 4 6 10 12],[2 4 6 10 12]);
Bmat_gir = Bmat([2 4 6 10 12],[1]);

Cmat_gir = [0 0 0 0 1];
Dmat_gir = 0;

[num_pade,den_pade] = pade(2.0,4);

for iiVs = 1:nVs
    
    Vs = Vsvect(iiVs)
	
    % Point de linearisation du modele
    Vs = Vsvect(iiVs)
	
	for jj=1:ALL.nb_states
		state = ALL.states{jj};
		Xeq.(state) = 0;
	end
 
	for jj=1:ALL.nb_inputs
		inputs = ALL.inputs{jj};
		XUeq.(inputs) = 0;
	end 
	
	Xeq.U      = Vs        ;
	Xeq.Cx0    = param.Cx0 ;
	XUeq.NT1   = CalculerEquilibreLOC(Xeq,XUeq,param,30.0) ;
	
	NT1 = XUeq.NT1;
	
	[Amat_copy, Bmat_copy] = subs_eq_gir(Amat_gir,Bmat_gir,s,u,seq,ueq, Vs, NT1);

	[Alat,Blat,Clat,Dlat] = CalculerModeleLinearise(Xeq, XUeq, param, GIR) ;
	
	Amat1 = eval(vpa(Amat_copy));
	Bmat1 = eval(vpa(Bmat_copy));

	Amat2 = Alat;
	Bmat2 = Blat;

	Amat1
	Amat2
	Amat1 - Amat2
	
end