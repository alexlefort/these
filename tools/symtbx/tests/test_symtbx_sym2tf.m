function bool = test_sym2tf

	syms s;
	
	f = (3*s^2 + 2*s + 1)/(4*s^6 + 3*s^3 + s^2 + 1);
	
	res      = symtbx_sym2tf(f,s);
	res_test = tf([3 2 1],[4 0 0 3 1 0 1]); 
	
	bool = isequaln(res_test,res);