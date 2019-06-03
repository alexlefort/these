function save_poly(c, z, k, p, name)

	exp_max = 12;
	
	kv  = fields(k);
	nkv = size(kv,1);
	
	expr = char(c);

	for ii=1:nkv
		ii
		kv{ii}
		char_x = strcat('x[',num2str(ii-1),']');
		expr   = strrep(expr ,kv{ii}, char_x);
	end

	pv  = fields(p);
	npv = size(pv,1);
	
	for ii=1:npv
		ii
		pv{ii}
		char_p = strcat('p[',num2str(ii-1),']');
		expr   = strrep(expr ,pv{ii}, char_p);
	end

	for ii=exp_max:(-1):1
		if ii==1
			char_w = 'w';
		else
			char_w = strcat('w^',num2str(ii));
		end
		char_p = strcat('exp(ln(10)*',num2str(ii),'*p[',num2str(npv),'])');
		expr = strrep(expr, char_w, char_p);
	end
	
	symvar(c)

	file = fopen(name,'w');
	
	mes = strcat('function f(x[' , num2str(nkv) , '],p[' , num2str(npv+1) , ']) \n return(' , expr , ') \n end');
	
	fprintf(file,mes);
	fclose(file);
