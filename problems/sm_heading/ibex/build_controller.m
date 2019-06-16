%%
%% Build a depth controller using stern plane with an integrator on Z :
%%
%% X = [Psi V R]
%% Y = [Psi V R]
%% U = [Alpha]

function model_ctrl = build_controller(Vs,g)

    %% Controller gains 
    if (nargin == 1)
        g.kpsi = sym('kpsi', 'real'); 
        g.kr   = sym('kr'  , 'real'); 
        g.ir   = sym('ir'  , 'real'); 
    end
    
    sat = sym('69813/1000000','real'); % Saturation en vitesse modelise comme un 1er ordre
    syms s;
    %% Transfert matrix representation
    sym_tf = [((g.kpsi*g.kr + g.ir)+g.ir*g.kpsi/s) 0 g.kr]*1/(s/sat + 1);
    
    %% State space representation     
    sym_ss.a = [0 0 ; sat*g.ir*g.kpsi -sat];
    sym_ss.b = [1 0 0 ; sat*(g.kr*g.kpsi +g.ir) 0 sat*g.kr];
    sym_ss.c = [0 1];
    sym_ss.d = [0 0 0];

    %% Save struct
    model_ctrl.gains = g;
    model_ctrl.symtf = sym_tf;
    model_ctrl.symss = sym_ss;

    model_ctrl.gnames.kpsi = sym('kpsi', 'real'); 
    model_ctrl.gnames.kr   = sym('kr'  , 'real'); 
    model_ctrl.gnames.ir   = sym('ir'  , 'real'); 