%% Horner Algorithm for multivariate polynomial : take
%%

function p_new = symtbx_horner(p,aux)

    if strcmp(aux,'simple')
        p_new = symtbx_horner_recurse(p);
    elseif strcmp(aux,'nodegree')
        pstr   = symtbx_deg2str(p);
        pstr_h = symtbx_horner_recurse(pstr);
        symbols = symvar(pstr_h);
        n = length(symbols);
        p_new = symtbx_str2deg(pstr_h);
    end
