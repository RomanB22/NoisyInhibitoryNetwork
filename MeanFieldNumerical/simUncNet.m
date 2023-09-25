%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    simulations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r0,x0]=simUncNet(ge0,sig0,Vre,Vth,EL,a,b,tau)

T=4000; 
dt=0.01; 
ns=ceil(T/dt); 
t=dt*[1:ns]';

V=zeros(ns,1);
x=zeros(ns,1);
r=zeros(ns,1);

V(1)=EL;    %initial conditions
x(1)=0;

for k=1:ns-1
 
    [I,IF]=currents(V(k),x(k),ge0,sig0,dt,tau);
    [xiV,txV]=activ_var(V(k),a,b);
    
    V(k+1)=V(k) + dt*(IF-I);  
    x(k+1)=x(k) + dt*(xiV-x(k))./txV;
        
    if V(k+1)>Vth
        r(k+1)=1;
        V(k+1)=Vre;
    end

end

trelax=1000;             % time to allow transients to decay
srelax=ceil(trelax/dt);

x0=mean(x(srelax:ns));   % calculate means with transients ignored
r0=sum(r(srelax:ns))/(t(ns)-trelax);
end