%%
%% Build a depth controller using stern plane with an integrator on Z :
%%
%% X = [Z Theta dZ dTheta]
%% Y = [Z Theta dZ dTheta]
%% U = [Beta1]

function model_ctrl = build_controller_simple(gains)

    %% Controller gains 
    if (nargin == 0)
        gains.kpt1  = sym('kpt1' ,'real'); 
        gains.kdt1  = sym('kdt1' ,'real');
        gains.kiz2  = sym('kiz2' ,'real');
        gains.kpz2  = sym('kpz2' ,'real'); 
        gains.kdz2  = sym('kdz2' ,'real'); 
    end
    
    syms s;
    
    %% Transfert matrix representation
    sym_tf = [0                            gains.kpt1 0          gains.kdt1  ; ...
              (gains.kpz2  + 0.1*gains.kiz2/s) 0          gains.kdz2 0          ];
    
    %% State space representation     
    sym_ss.a = zeros(1);
    sym_ss.b = [1 0 0 0];
    sym_ss.c = [0 ; 0.1*gains.kiz2];
    sym_ss.d = [0 gains.kpt1 0 gains.kdt1 ; ...
                gains.kpz2 0 gains.kdz2 0];

    %% Save struct
    model_ctrl.gains = gains;
    model_ctrl.symtf = sym_tf;
    model_ctrl.symss = sym_ss;

    model_ctrl.gnames.kpt1  = sym('kpt1' ,'real'); 
    model_ctrl.gnames.kdt1  = sym('kdt1' ,'real');
    model_ctrl.gnames.kiz2  = sym('kiz2' ,'real');
    model_ctrl.gnames.kpz2  = sym('kpz2' ,'real'); 
    model_ctrl.gnames.kdz2  = sym('kdz2' ,'real'); 

	