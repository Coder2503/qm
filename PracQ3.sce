clc;clear;clf();
h=197.3;m=940;k=100;
rmin=1e-10;rmax=10
n=1000;b=[0;10;30]
r=linspace(rmin,rmax,n);
d=r(2)-r(1);
t=-(h^2)/(2*m*(d^2));
V=zeros(n,n);
U=-2*eye(n,n)+diag(ones(n-1,1),1)+diag(ones(n-1,1),-1);
for j=1:3
    for i=1:n
        V(i,i)=0.5*k*(r(i)^2)+(b(j)*(r(i)^3))/3;
    end
    H=t*U+V;
    [U1,EV]=spec(H);
    E=round(diag(EV)*1000)/1000;
    disp("Value of b : "+string(b(j)))
    disp("Ground state energy : "+string(E(1))+" MeV")
    figure(j-1);scf(j-1);clf(j-1)
    plot(r',U1(:,1),'linewidth',2)
    title("Plot of wavefunction for V=0.5*k*r^2 +("+string(b(j))+"*r^3)/3")
    xlabel("r(A)")
    ylabel("Wavefunction")
    a=gca()
    a.x_location="origin";
    a.y_location="origin";
    xgrid()
    legend("Ground State","First excited state",4)
end
