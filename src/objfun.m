function S=objfun(p,X,Y,Z)

S=sum(sum((p(1)*exp(-0.5/p(4)^2*(cos(p(6))*(X-p(2))+sin(p(6))*(Y-p(3))).^2-0.5/p(5)^2*(-sin(p(6))*(X-p(2))+cos(p(6))*(Y-p(3))).^2)+p(7)-Z).^2));

