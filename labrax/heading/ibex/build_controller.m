%%
%% Build a depth controller using stern plane with an integrator on Z :
%%
%% X = [Z Theta dZ dTheta]
%% Y = [Z Theta dZ dTheta]
%% U = [Beta1]

function model_ctrl = build_controller(Vs,g)

    %% Controller gains 
    if (nargin == 1)
        g.kz  = sym('kz' , 'real'); 
        g.kpi = sym('kpi', 'real'); 
        g.ipi = sym('ipi', 'real'); 
        g.kq  = sym('kq' , 'real');
    end
    
    syms s;
    r = (g.kpi + g.ipi/s);
    %% Transfert matrix representation
    sym_tf = [g.kz*r r -r/Vs g.kq];
    
    %% State space representation     
    sym_ss.a = zeros(1);
    sym_ss.b = [g.kz*g.ipi g.ipi -g.ipi/Vs g.ipi];
    sym_ss.c = 1;
    sym_ss.d = [g.kz*g.kpi g.kpi -g.kpi/Vs g.kq];

    %% Save struct
    model_ctrl.gains = g;
    model_ctrl.symtf = sym_tf;
    model_ctrl.symss = sym_ss;

    model_ctrl.gnames.kz   = sym('kz' , 'real'); 
    model_ctrl.gnames.kpi  = sym('kpi', 'real'); 
    model_ctrl.gnames.ipi  = sym('ipi', 'real'); 
    model_ctrl.gnames.kq   = sym('kq' , 'real');
	