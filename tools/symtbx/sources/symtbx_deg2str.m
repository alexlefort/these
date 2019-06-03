%% For the multivariate polynomial p, replace each monomial as :
%% x^6 -> x6

function pstr = symtbx_deg2str(p)

    deg_max = symtbx_degmax(p);
    symbols = symvar(p);
    n = length(symbols);
    pstr = expand(p);

    for ii=deg_max:-1:1 %% Reverse loop in case of degrees > 10
        for jj=1:n
            nstr = strcat(char(symbols(jj)),num2str(ii));
            str = sym(nstr);
            pstr = subs(pstr,symbols(jj)^ii,str);
        end
    end