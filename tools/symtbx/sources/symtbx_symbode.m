%% return bode module from a symbolic tf

function res = symtbx_symbode(trans, s, w)

    var = symvar(trans);

    if (length(var) > 1)
        error('symtbx_symbode : undefined parameters');
    end

    res = 20*log10(abs(eval(subs(trans,s,w))));
    
