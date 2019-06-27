%% Take a structure of vector values and build a table with 
%% all possible combination

function res = param_org_table(p)

f = fields(p);
nparam = length(f);

for ii=1:nparam
    c{ii} = p.(f{ii});
end
b = cell(1,numel(c));
[b{:}] = ndgrid(c{:});

for ii=1:nparam
    res(:,ii) = reshape(b{ii},1,[]);
end
[n,~] = size(res);

disp([' param_org_table : number of points : ' num2str(n)]);