clear all

m = 0.0123; %rapporto masse corpo celeste 1 / corpo celeste 2
m2 = m/(1+m); %le masse sono normalizzate in modo che m1+m2=1
m1 = 1-m2;
D = 3.844e5; %distanza tra i corpi celesti in km


%dati iniziali
distorbit = 59699; %in km, distanza dal centro del corpo 1 dell'astronave
x0 = -m2+distorbit/D*(m2+m1); %normalizzazione
y0 = 0; 
u0 = 0;

d1 = sqrt((x0+m2)^2+y0^2); %distanza corpo1-astronave
d2 = sqrt((x0-m1)^2+y0^2); %distanza corpo2-astronave


%inizializziamo il problema con integrale di Jacobi J
J0 = -3.17948;
v0 = sqrt(J0-u0^2+x0^2+y0^2+2*(m1/d1+m2/d2));
z0 = [x0,y0,u0,v0]; %punto iniziale

%calcolo dell'orbita con metodo odesolver
opts = odeset('AbsTol',eps,'RelTol',2.5e-14);
endtime = 7000; %tempo finale a cui arrestare le computazioni
tspan = [0,endtime];
tempi = readmatrix('TempiPoincare.txt');
sol = table2array(readtable("SezionePoincare.txt"))';
scatter(sol(1,:),sol(2,:),'.','k')
a = [x0,0.2];
b = [1,0];
v = (sol(1,:)-a(1)).^2+(sol(2,:)-a(2)).^2;
[l,index_a] = min(v);
index_a = min(index_a);
v = (sol(1,:)-b(1)).^2+(sol(2,:)-b(2)).^2;
[l,index_b] = min(v);
index_b = max(index_b);
[index_a,index_b]
del = 0.03;
epp = 0.01;
loop_length = 12;
%[final,final_index] = min(abs(sol.x-sol.x(find(sol.y==b))-loop_length));

array = sol(:,min(index_a,index_b):max(index_a,index_b));

[f,t] = BolltMeissAlgo(array,del,epp,loop_length,index_a,index_b,tempi,sol);

hold on
plot(-m2,0,['x','b']); %in blu il corpo 1
plot(m1,0,['x','r']); %in rosso il corpo 2
scatter(f(1,:),f(2,:),'.')
for o=1:length(f)
    indicideitempi(o)=indice(sol,f(:,o));
end





