m = 0.0123; %rapporto masse corpo celeste 1 / corpo celeste 2
m2 = m/(1+m); %le masse sono normalizzate in modo che m1+m2=1
m1 = 1-m2;
D = 3.844e5; %distanza tra i corpi celesti in km


%dati iniziali
%x0 = 0.26
distorbit = 59669; %in km, distanza dal centro del corpo 1 dell'astronave
x0 = -m2+distorbit/D*(m2+m1); %normalizzazione
y0 = 0; 
u0 = 0;

d1 = sqrt((x0+m2)^2+y0^2); %distanza corpo1-astronave
d2 = sqrt((x0-m1)^2+y0^2); %distanza corpo2-astronave



%inizializziamo il problema con integrale di Jacobi J0
J0 = -3.17948; %integrale di Jacobi dell'orbita (costante del moto)
v0 = sqrt(J0-u0^2+x0^2+y0^2+2*(m1/d1+m2/d2));
z0 = [x0,y0,u0,v0]; %punto iniziale

%iniziamo a mappare sul grafico la posizione dei due corpi celesti
hold on
plot(-m2,0,['x','b']); %in blu il corpo 1
plot(m1,0,['x','r']); %in rosso il corpo 2
%calcolo dell'orbita con metodo odesolver
opts = odeset('AbsTol',eps,'RelTol',2.5e-14);
endtime = 300000; %tempo finale a cui arrestare le computazioni
nsteps = 1000; %numero di passi in cui spezzare il vettore dei tempi (metodo per monitorare il progresso delle computazioni)
nsteps = fix(max(1,nsteps)); %in questo modo si evitano error
t = 0;
k = 0;
for j = 1:nsteps
    fprintf("\n%2.1f%%",(j-1)/nsteps*100); % stampa percentuale di progresso

    %calcolo della soluzione ai passi successivi, si consiglia ode113 o ode45
    sol = ode113(@(t,y)ThreeBodyRes(t,y,m),[(j-1)*endtime/nsteps,j*endtime/nsteps],z0,opts);
    
    %calcolo del dato iniziale del prossimo intervallo di tempo parziale
    l = length(sol.y); %lunghezza del segmento di soluzione attuale
    z0 = sol.y(:,l); %dato iniziale successivo

    %costruzione del grafico sul piano di Poincaré y=0,v>0 e check
    %della conservazione dell'integrale di Jacobi
    tol = 1e-6; %tolleranza per l'errore relativo per integrale Jacobi
    for i=2:l
        %se il punto (x,u) appartiene al piano di Poincaré viene salvato in f
        if (sol.y(2,i-1)*sol.y(2,i)<0)&&(sol.y(4,i)>0)&&(abs(sol.y(1,i))<10)
            k = k+1;
            f(1,k)=sol.y(1,i);
            f(2,k)=sol.y(3,i);
            g(k)=sol.x(i);
        end

        %calcolo dell'integrale di Jacobi al punto i-simo dell'intervallo
        %temporale attuale
        d1 = sqrt((sol.y(1,i)+m2)^2+sol.y(2,i)^2);
        d2 = sqrt((sol.y(1,i)-m1)^2+sol.y(2,i)^2);
        J = sol.y(3,i)^2+sol.y(4,i)^2-sol.y(1,i)^2-sol.y(2,i)^2-2*(m1/d1+m2/d2);

        %check dell'errore relativo dell'integrale di Jacobi rispetto a
        %quello iniziale.
        if abs((J-J0)/J0)>tol
            fprintf("\nATTENZIONE! \nL'integrale di Jacobi si conserva con errore relativo minore di %d solo per i primi %d passi",[tol,i]);
            fprintf("\nIl totale dei passi è %d",l);
            break
        end
    end  
    
    %grafico delle intersezioni col piano di Poincaré
    %scatter(f(1,:),f(2,:),['.','k']);
    %scatter(sol.y(1,:),sol.y(2,:),['.','k']);
end

fprintf("\n100%%");

x = f(1,:)';
u = f(2,:)';

T = table(x,u)
writetable(T,'SezionePoincare.txt');

writematrix(g,'TempiPoincare.txt')














