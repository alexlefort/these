%% For a multivariate polynomial, find the higher degree of monomial
%%

function deg = symtbx_degmax(p)

    symbols = symvar(p);
    n = length(symbols);

    deg = 0;
    for ii = 1:n
        deg = max(deg, symtbx_get_sympoly_degree(p,symbols(ii)));
    end