function [M]=gm(l,r)

% [M]=gm(l,r). M is a matrix of geometric moments up to the r-th order of the image l
% l(n1,n2) - image matrix
% \mu{pq} = M(p+1,q+1)

 [n1,n2]=size(l);

 a=linspace(1,n2,n2)-1;
 c=linspace(1,n1,n1)-1;

for i=1:r+1;
for j=1:r+1;
p=i-1;
q=j-1;
A=a.^p;
C=c.^q;
M(i,j)=C*l*A';
end;
end;

