function i = indice(sol,val)
v = (sol(1,:)-val(1)).^2+(sol(2,:)-val(2)).^2;
[l,i] = min(v);
i = min(i);
end
