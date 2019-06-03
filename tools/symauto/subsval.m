function  [res] = subsval(inp, p, psub)

	pf = fields(psub);
	
	nb_p = size(pf,1);

    res = inp;
	for ii=1:nb_p
		res = subs(res,p.(pf{ii}) , psub.(pf{ii}));
	end