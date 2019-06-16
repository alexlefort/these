%%
%% Build a depth controller using stern plane with an integrator on Z :
%%
%% X = [Z Theta dZ dTheta]
%% Y = [Z Theta dZ dTheta]
%% U = [Beta1]

function model_filter = build_filter(g)

    %% Controller gains 
    if (nargin == 0)
        g.oz  = sym('oz','real'); 
        g.ot  = sym('ot','real'); 
        g.kz  = sym('kz','real'); 
        g.kt  = sym('kt','real');
    end
    
    syms s;
    

    %% Filter
    Fz = 1/(s^2/g.oz^2 + 2*g.kz/g.oz +1);
    Ft = 1/(s^2/g.ot^2 + 2*g.kt/g.ot +1);
       
    %% Transfert matrix representation
    sym_tf = [Fz    0   ; ...
               0    Ft  ; ...
              Fz*s  0   ; ...
               0    Ft*s];

           
    A = [      0            1       0            0 ; ...
         -g.oz^2 -2*g.kz*g.oz       0            0 ; ...
               0            0       0            1 ; ...
               0            0 -g.ot^2 -2*g.kt*g.ot];

    B = [     0          0 ; ...
         g.oz^2          0 ; ...
              0          0 ; ...
              0    g.ot^2 ];
          
    C = eye(4);
    D = zeros(4,2);
    
    %% State space representation     
    sym_ss.a = A;
    sym_ss.b = B;
    sym_ss.c = C;
    sym_ss.d = D;

    %% Save struct
    model_filter.gains = g;
    model_filter.symtf = sym_tf;
    model_filter.symss = sym_ss;

    model_filter.gnames.oz    = sym('oz','real'); 
    model_filter.gnames.ot    = sym('ot','real'); 
    model_filter.gnames.kz    = sym('kz','real'); 
    model_filter.gnames.kt    = sym('kt','real');
