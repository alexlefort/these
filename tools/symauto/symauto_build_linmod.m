function [Amat, Bmat] = symauto_build_linmod(dx,s,u,seq,ueq)
	
	ns = fields(s);
	nu = fields(u);
	
	nseq = fields(seq);
	nueq = fields(ueq);
	
	nb_s = size(ns,1);
	nb_u = size(nu,1);
	
	Amat = sym('Amat',[nb_s nb_s]);
	Bmat = sym('Bmat',[nb_s nb_u]);
	
	
	%% Define eq states and input

	for ii=1:nb_s
		for jj=1:nb_s
			aux1 = s.(ns{ii});
			aux2 = s.(ns{ii}) + seq.(nseq{ii});
			dx.(ns{jj}) = subs(dx.(ns{jj}),{aux1}, {aux2});
		end
	end
	
	for ii=1:nb_u
		for jj=1:nb_s
			aux1 = u.(nu{ii});
			aux2 = u.(nu{ii}) + ueq.(nueq{ii});
		
			dx.(ns{jj}) = subs(dx.(ns{jj}),{aux1}, {aux2});
		end
	end

	%% Compute the Jacobian
	
	for ii=1:nb_s
		for jj=1:nb_s
			aux1 = s.(ns{jj});
			Amat(ii,jj) = diff(dx.(ns{ii}),aux1);
		end
	end
	
	
	for ii=1:nb_s
		for jj=1:nb_u
			aux1 = u.(nu{jj});
			Bmat(ii,jj) = diff(dx.(ns{ii}),aux1);
		end
    end
	
	%% Erease non steady states 
	
	for ii=1:nb_s
		Amat = subs(Amat,{s.(ns{ii})},{0});
		Bmat = subs(Bmat,{s.(ns{ii})},{0});
	end
	
	
	for ii=1:nb_u
		Amat = subs(Amat,{u.(nu{ii})},{0});
		Bmat = subs(Bmat,{u.(nu{ii})},{0});
	end
