%% For the multivariate polynomial p, replace each monomial as :
%% x6 -> x^6

function pdeg = symtbx_str2deg(p)

    symbols = symvar(p);
    n = length(symbols);
    pdeg = p;

    for ii=1:n
        strr   = char(symbols(ii));
        name   = strr(1:end-1);
        order  = strr(end);
        symbol = sym(name);
        pdeg = subs(pdeg,strr,symbol^order);
    end