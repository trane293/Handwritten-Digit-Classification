function [A]=zm(img,rd,norm)
% [A]=zm(img,rd,norm) computes Zernike moments up the rd-th order of the image img
% img(n1,n2) - image matrix
% if norm=0, no normalization to scaling will be provided
% if norm=1, moments will be normalized to scaling (but not to density of sampling), 
% if norm=2, origin of coordinates is in the center of the image, else in
% the centroid
% The moment A_{mn} = A(m+1,(m-n)/2+1)

 [n1,n2]=size(img);

m00 =sum(sum(img));
if norm==2
    tx=(n2-1)/2;
    ty=(n1-1)/2;
else
	w=linspace(1,n2,n2);
	v=linspace(1,n1,n1);
	if m00~=0
        tx=(sum(img*w'))/m00; 
        ty=(sum(v*img))/m00;
	else
        tx=0;
        ty=0;
	end
end
if norm==0
    rmax=max([(1-tx)^2+(1-ty)^2,(n2-tx)^2+(1-ty)^2,(n2-tx)^2+(n1-ty)^2,(1-tx)^2+(n1-ty)^2]);
    rmax=sqrt(rmax);
elseif norm==2
    rmax=sqrt(n1*n2)*sqrt(n2/n1+n1/n2)/2;
else
    rmax=sqrt(m00)*sqrt(n2/n1+n1/n2)/2;
end

[x,y] = meshgrid(1:n2,1:n1);
x=(x-tx)/rmax;
y=(y-ty)/rmax;
r=sqrt(x.^2+y.^2);
theta=atan2(y,x);

%Kintner method
for n=-rd:rd  %repetition
    an=abs(n);
    rmn0=r.^an;
    vmn=rmn0.*exp(i*n*theta);
    prodiv=img.*conj(vmn);
    A(an+1,(an-n)/2+1)=(an+1)/pi*sum(sum(prodiv));
    if(rd-an>=2)
        rmn2=(an+2)*r.^(an+2)-(an+1)*r.^an;
        vmn=rmn2.*exp(i*n*theta);
        prodiv=img.*conj(vmn);
        A(an+3,(an+2-n)/2+1)=(an+3)/pi*sum(sum(prodiv));
    end
    for m=an+4:2:rd %order
        k1=(m+n)*(m-n)*(m-2)/2;
        k2=2*m*(m-1)*(m-2);
        k3=-n^2*(m-1)-m*(m-1)*(m-2);
        k4=-m*(m+n-2)*(m-n-2)/2;
        rmn4=((k2*r.^2+k3*ones(n1,n2)).*rmn2+k4*rmn0)/k1;
        vmn=rmn4.*exp(i*n*theta);
        prodiv=img.*conj(vmn);
        A(m+1,(m-n)/2+1)=(m+1)/pi*sum(sum(prodiv));
        rmn0=rmn2;
        rmn2=rmn4;
    end
end