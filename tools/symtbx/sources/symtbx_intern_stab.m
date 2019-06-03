%% return internal stability polynomial of a closed-loop system

function [poly, mat] = symtbx_intern_stab(sys, sys_ctrl)
    
    A = sys.a;
    B = sys.b;
    C = sys.c;
    %%D = sys.D;

    Ak = sys_ctrl.a;
    Bk = sys_ctrl.b;
    Ck = sys_ctrl.c;
    Dk = sys_ctrl.d;

    if (Ak == 0 && Ck == 0)
        mat = (A+B*(-Dk)*C);
    else
        mat = [(A+B*(-Dk)*C) (B*Ck) ; ...
               (-Bk*C)       (Ak)  ];
    end

    syms s;

    [n,m] = size(mat);
    poly  = simplify(det(s*eye(n) - mat));
