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
    
    syms s;
    %% Transfert matrix representation
    sym_tf = [((g.kpsi*g.kr + g.ir)+g.ir*g.kpsi/s) 0 g.kr];
    
    %% State space representation     
    sym_ss.a = zeros(1);
    sym_ss.b = [1 0 0];
    sym_ss.c = g.ir*g.kpsi;
    sym_ss.d = [(g.kr*g.kpsi +g.ir) 0 g.kr];

    %% Save struct
    model_ctrl.gains = g;
    model_ctrl.symtf = sym_tf;
    model_ctrl.symss = sym_ss;

    model_ctrl.gnames.kpsi = sym('kpsi', 'real'); 
    model_ctrl.gnames.kr   = sym('kr'  , 'real'); 
    model_ctrl.gnames.ir   = sym('ir'  , 'real'); 