clear all

p.p1 = -0.235:0.02:0.235;
p.p2 = 0.998:0.0005:1.002;
p.p3 = 5:1:10;

f = fields(p);
nparam = length(f);

for ii=1:nparam
    c{ii} = p.(f{ii});
end
b = cell(1,numel(c));
[b{:}] = ndgrid(c{:});
disp(b);

for ii=1:nparam
    res(:,ii) = reshape(b{ii},1,[]);
end

disp(size(res));