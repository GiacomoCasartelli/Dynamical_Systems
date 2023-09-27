function [f,times] = BolltMeissAlgo(array,del,epp,loop_length,index_a,index_b,tempi,sol)
 
opts = odeset('AbsTol',eps,'RelTol',2.5e-14);
mass = 0.0123;
m2 = mass/(1+mass); 
m1 = 1-m2;
l = length(sol);
check = false;
k2 = length(array);
Booleanoss=true;
for k=1+loop_length:k2-loop_length 
    %fprintf("\nStart of loop number %d of %d",k-k1+1-loop_length,k2-k1-2*loop_length+1)
    for j=k2-loop_length:-1:k+loop_length %we want to find zj with del-recurrence
        %fprintf("\n%2.2f%%",(1-(j-k)/(l-k))*100)
        if (norm(array(:,j)-array(:,k))<del) %check for recurrence
            %we search for a patch
            %let's build a new orbit
            if (k>k2/2)&&Booleanoss
                fprintf("\nwe are making progress%d",j)
                Booleanoss=false;
            end
            zi = array(:,k-loop_length);
            zj = array(:,j+loop_length);
            fu = [0;1];
            t = newton(@(s)dot((sol(:,min(l,indice(sol,zi+s.*fu)+2*loop_length))- zj),fu),.5,1,1e-5,100);
            z0 = zi + t.*fu;
            index = indice(sol,z0);
            index_1 = max(1,index);
            index_2 = min(l,index+2*loop_length);
            neworbit_1 = sol(:,index_1);
            neworbit_2 = sol(:,index_2);
            %check if it is a valid patch    
            if (norm(zj-neworbit_2)<epp)&&(norm(zi-neworbit_1)<epp)
                    fprintf("\nA suitable patch was found")
                    check = true;
                    times = [k-loop_length,j+loop_length,index_1,index_2]
                    ind_3 = k-loop_length;
                    f = zeros(2,k2-(j-k)-2*loop_length+(index_2-index_1)+1);
                    f(:,1:ind_3) = array(:,1:ind_3);
                    f(:,ind_3+1:ind_3+1+index_2-index_1) = sol(:,index_1:index_2);
                    f(:,ind_3+2+index_2-index_1:k2-(j-k)-2*loop_length+(index_2-index_1)+1) = array(:,j+loop_length+1:k2);
                    tempo_per_patch = tempi(index_2)-tempi(index_1);
                    fprintf("\nThe patch time is %f",tempo_per_patch)
                    fprintf("\nThe patch times are %f and %f",tempi(index_1),tempi(index_2))
                    break
            end
        end
        if check
            break
        end
    end
    if check
        break
    end
end

if check
    [f,times] = BolltMeissAlgo(f,del,epp,loop_length,index_a,index_b,tempi,sol);
else
    f=array;
    times=[index_a,index_b];
end

end

