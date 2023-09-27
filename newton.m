function g = newton(f,x0,x1,tol,Nmax)

if nargin < 3
    error("\nNot enough inputs.")
end
if nargin==3
    tol = 1e-6;
    Nmax = 10000;
    fprintf("\ntol set to %f",tol)
    fprintf("\nNmax set to %f",Nmax)
end
if nargin==4
    Nmax = 10000;
    fprintf("\nNmax set to %f",Nmax)
end

g1 = x1;
g0 = x0;
for n=2:Nmax
    g = g1 - f(g1).*(g1-g0)./(f(g1)-f(g0));   
    if abs(g-g1)<tol
        break
    end
    g0 = g1;
    g1 = g;
end


