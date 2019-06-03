function  res = symtbx_subs_struct(x, p, psub)

	pf = fields(psub);
	nb_p = size(pf,1);

	res = x;

	for ii=1:nb_p
		res = subs(res, p.(pf{ii}), psub.(pf{ii}));
	end