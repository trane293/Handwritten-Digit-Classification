function [A]=zmgm(gm,r)
% [A]=zmgm(gm,r)
% Conversion of geometric moments gm to Zernike moments A up to the r-th order 
% A_{nm} is in A(n+1,(m+n)/2+1)

%Prata method
rmn=zeros(r+1,r+1,r+1);
for m=0:r
    rmn(m+1,1,1)=1;
    rmn(m+1,m+1,1)=1;
end
for m=1:r
    for n=2-m:2:m-2
        rmn(m+1,(m-n)/2+1,1:r+1)=2*m/(m+n)*rmn(m,(m-n)/2+1,1:r+1);
        rmn(m+1,(m-n)/2+1,3:r+1)=rmn(m+1,(m-n)/2+1,3:r+1)-(m-n)/(m+n)*rmn(m-1,(m-2-n)/2+1,1:r-1);
    end
end

A=zeros(r+1)+i*zeros(r+1);
for n=0:r;
    for l=-n:2:n;
        if(l>0)
            w=-i;
        else
            w=i;	%imaginary unit
        end
        la=abs(l);
        for k=la:2:n;
            q=(k-la)/2;
            for j=0:q;
                for m=0:la;
                    A(n+1,(l+n)/2+1)=A(n+1,(l+n)/2+1)+(n+1)/pi*w^m*...
                        nchoosek(q,j)*nchoosek(la,m)*rmn(n+1,(n-l)/2+1,n-k+1)*gm(k-2*j-m+1,2*j+m+1);
                end
            end
        end
    end
end

