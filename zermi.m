function [invts,invmat,ind,normom]=zermi(a,r,thres,typecomp,typeinv)
% [invts,invmat,ind,normom]=zermi(a,r,thres,typecomp,typeinv) calculates 
% values of rotation Zernike invariants of the image a (no im[l. value) 
% up to the order r (impl. 3), 
% thres is threshold of non-zero moments, impl. 1e-3
% typecomp is method of computation: if typecomp=0, moments are
% computed directly by Kintner method (impl.), if typecomp=1, they are 
% converted from geometric ones,
% typeinv is type of rotation invariancy: if typeinv=0 (impl.), rotation 
% normalization is used, if typeinv=1, phase cancellation by multiplication 
% is used.
% invts  - values of the invariants
% invmat - values of the invariants in matrix form,
%            A_{r,ell}=invmat(r+1,(r+ell)/2+1);
% ind    - vector of indicators, ind(1,i)=i is the number of i-th invariant,
%          ind(2,i) is order of i-th invariant, ind(3,i) is repetition and
%          ind(4,i)=0 for real part and ind(4,i)=1 for imaginary part.
% normom - vector of three numbers describing the moment normalizing rotation:
%          order, repetition and normalizing angle.

if nargin<2
    r=3;
end
if nargin<3
    thres=1e-3;
end
if nargin<4
    typecomp=0;
end
if nargin<5
    typeinv=0;
end

[x1,y1]=find(a);
szx=size(x1);
szy=size(y1);
if(szx(1)>0 && szy(1)>0)
	mn1=min(x1);
	mx1=max(x1);
	mn2=min(y1);
	mx2=max(y1);
    a=a(mn1:mx1,mn2:mx2);
end
if typecomp==0
	am=zm(a,r,1);   %computation of Zernike moments
	am=am/am(1,1);  %normalization of moments to density of sampling + contrast
else
	[qm,pm]=meshgrid(0:r);
	gm=cm(a,r);                     %computation of central moments TRANSLATION INVARIANCE
	gm00=gm(1,1);
	gm=gm./(gm00.^((pm+qm+2)/2));   %normalization of moments to scaling SCALING INVARIANCE
	am=zmgm(gm,r);                  %conversion of moments to Zernike ones
end
mt=0;
nt=0;
theta=0;
flag=0;         %normalization of moments to rotation
for k1=1:r
    for k2=k1:2:r
        if k2>1
            if abs(am(k2+1,(k2+k1)/2+1))>thres
                theta=angle(am(k2+1,(k2+k1)/2+1))/k1;   %normalizing moment found
                flag=1;
                mt=k2;
                nt=k1;
            end
        end
        if flag
            break
        end
    end
    if flag
        break
    end
end
if flag
	for k1=2:r
        for k2=-k1:2:k1
            if typeinv==0
                am(k1+1,(k1+k2)/2+1)=am(k1+1,(k1+k2)/2+1)*exp(-1i*k2*theta);
            else
                if k1~=mt || k2~=nt
                    am(k1+1,(k1+k2)/2+1)=am(k1+1,(k1+k2)/2+1)/am(mt+1,(mt+nt)/2+1)^(k2/nt);
                end
            end
        end
    end
    if typeinv~=0
        am(mt+1,(mt+nt)/2+1)=abs(am(mt+1,(mt+nt)/2+1));
    end
end
%end of normalization of moments to rotation

% invariants
pinv=1;
for k1=2:r
    for k2=k1:-2:0
        invts(pinv)=real(am(k1+1,(k1+k2)/2+1));
        cmp(pinv,:)=[pinv,k1,k2,0];
        pinv=pinv+1;
        if k2~=0 && (k1~=3 || k2~=1)
            invts(pinv)=imag(am(k1+1,(k1+k2)/2+1));
            cmp(pinv,:)=[pinv,k1,k2,1];
            pinv=pinv+1;
        end
    end
end
invmat=am;
ind=cmp';
normom=[mt,nt,theta];
