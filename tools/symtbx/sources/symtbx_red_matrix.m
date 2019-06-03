%% Give reduced matrix 
%%

function M_red = symtbx_red_matrix(M,row,col)

	[nrows,ncols] = size(M);

	if (row>nrows && col>ncols)
		error('invalid rol or col index i red_matrix');
		M_red = 0;
	else
		M_red = M([1:row-1 row+1:end], [1:col-1 col+1:end]);
	end