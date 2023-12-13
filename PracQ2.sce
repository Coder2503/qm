clear;clc;
h=1973;m=0.511e6;e=3.795
rmin=1e-10;rmax=10;
n=1000;
a0=[3;5;7]
r=linspace(rmin,rmax,n);
d=r(2)-r(1);
t=-(h^2)/(2*m*(d^2));
V=zeros(n,n);
K=-2*eye(n,n)+diag(ones(n-1,1),1) +diag(ones(n-1,1),-1);
for k=1:3
for i=1:n
    V(i,i)=-(e^2)*%e^(-r(i)/a0(k))/r(i);
end
H=t*K+V;
[U1,EV]=spec(H);
E=round(diag(EV)*1000)/1000;
disp("Value of a = "+string(a0(k)))
disp("Ground state energy : "+string(E(2))+" eV")
figure(k);scf(k);clf(k);
plot(r',U1(:,2),'linewidth',2)
legend('Ground state',4)
title("Plot of wavefunction for e^(2)/r *e^(-r/"+string(a0(k))+") Potential",'fontsize',5)
xlabel("r(A)",'fontsize',4)
ylabel("Wavefunction",'fontsize',4)
xgrid()
a=gca();
a.x_location="origin";
a.y_location="origin";
end
