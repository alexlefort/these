%% Return the transfert matrix from a state space representation

function sym_tf = symtbx_ss2symtf(ss_t)

    syms s;

    A = ss_t.a;
    B = ss_t.b;
    C = ss_t.c;
    D = ss_t.d;

    [n,m] = size(A);
    sym_tf = C*inv(s*eye(n)-A)*B + D;