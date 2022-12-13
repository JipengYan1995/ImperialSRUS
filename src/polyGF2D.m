function [A,cx,cy,sx,sy,sita]=polyGF2D(x,y,z,h)

f=[1 2 1;2 4 2;1 2 1]/16; 
z=filter2(f,z);    
for k=1:1    
    z=filter2(f,z);
end
zmax=max(max(z));        
xnew=x(find(z>(zmax*h)));
ynew=y(find(z>(zmax*h)));
znew=z(find(z>(zmax*h)));
zlognew=log(znew);
a = zeros(max(size(xnew)),6);
n=1;
for i1= 0:2
   for j1=0:2-i1
       a(:,n) = (xnew.^i1).*(ynew.^j1);
       n=n+1;
   end
end
p = (a\zlognew);
c(1)=-p(6);c(2)=-p(3);c(3)=-p(5);c(4)=-p(4);c(5)=-p(2); c(6)=p(1);
sitap=0.5*acot((c(1)-c(2))/c(3));
sxp=sqrt(1/((c(1)-c(2))/cos(2*sitap)+c(1)+c(2)));
syp=sqrt(1/(-(c(1)-c(2))/cos(2*sitap)+c(1)+c(2)));

MA=([-2*(cos(sitap)^2/sxp^2+sin(sitap)^2/syp^2),sin(2*sitap)*(1/syp^2-1/sxp^2);...
    sin(2*sitap)*(1/syp^2-1/sxp^2),-2*(sin(sitap)^2/sxp^2+cos(sitap)^2/syp^2)]);
MB=2*[c(4),c(5)]';
cp=MA\MB;
cxp=cp(1);
cyp=cp(2);
Ap=exp(c(6)+(cos(sitap)*cxp+sin(sitap)*cyp).^2/2/sxp^2+...
  (-sin(sitap)*cxp+cos(sitap)*cyp).^2/2/syp^2);
A=Ap;cx=cxp;cy=cyp;sx=sxp;sy=syp;sita=sitap;