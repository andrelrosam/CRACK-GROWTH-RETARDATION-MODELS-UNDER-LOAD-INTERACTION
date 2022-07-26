function [s_i,itr] = gauss_seidel(S,f,L,alpha,s_i,s_f,g,n_zone)
tol = 0.02*s_f;
itr = 0;

n=length(s_i);
normVal=Inf; 

while normVal>tol
s_i_old=s_i;
    
    for i=1:n
        
        sigma=0;
        
        for j=1:i-1
                sigma=sigma+g(i,j)*s_i(j);
        end
        
        for j=i+1:n
                sigma=sigma+g(i,j)*s_i_old(j);
        end
        
        s_i(i)=(1/g(i,i))*(S*f(i)-L(i)-sigma);
        
        if i<=n_zone
            if s_i(i)>(alpha*s_f)
                s_i(i) = alpha*s_f;
            elseif s_i(i)<(-1*s_f)
                s_i(i) = (-1*s_f);
            end
        else
            if s_i(i)>0
                s_i(i) = 0;
            elseif s_i(i)<(-1*s_f)
                s_i(i) = (-1*s_f);
            end
        end
        
    end
    
    itr=itr+1;
    normVal=norm(s_i_old-s_i);
end