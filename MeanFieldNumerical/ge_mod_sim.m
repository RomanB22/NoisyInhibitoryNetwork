function [r0,ar,thr,ax,thx]=ge_mod_sim(fkHz,T,dt1,N,EL,ge0,ge1,Vre,Vth,sig0Sim,a,b,tau)

factor=1;
period=1/fkHz;                 % period in ms
nspp=ceil(period/dt1/factor);         % number of steps per period
dt=period/nspp;                % the actual dt
np=ceil(T/period);             % the number of periods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% oscillatory conductance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tp=dt*[0:nspp-1]';
ge1t=ge1*sin(2*pi*tp/period);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare the arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r=zeros(nspp,1);
x=zeros(nspp,1);

nE=10;                  % number of examples
VE=zeros(nspp,nE);
xE=zeros(nspp,nE);

Vold=EL*ones(1,N); 
[xold,dummy]=activ_var(Vold,a,b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a short equilibriation run without osc drive  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ns_therm=ceil(1000/dt);     

for k=1:ns_therm
    
    [I,IF]=currents(Vold,xold,ge0,sig0Sim,dt,tau);
    [xiV,txV]=activ_var(Vold,a,b);
                
    Vnew=Vold + dt*(IF-I);
    xnew=xold + dt*(xiV-xold)./txV;
                
    spiked=find(Vnew>Vth);              
    Vnew(spiked)=Vre;
   
    Vold=Vnew;            % moves the system on
    xold=xnew;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop over all the periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=1:np

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loop over a period
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for s=1:nspp

        [I,IF]=currents(Vold,xold,ge0+ge1t(s),sig0Sim,dt,tau);
        [xiV,txV]=activ_var(Vold,a,b);

        Vnew=Vold + dt*(IF-I);
        xnew=xold + dt*(xiV-xold)./txV;

        spiked=find(Vnew>Vth);                                  
        r(s)=r(s)+length(spiked);
        Vold(spiked)=Vth;
        Vnew(spiked)=Vre;

        VE(s,:)=Vold(1:nE);
        xE(s,:)=xold(1:nE);

        Vold=Vnew;
        xold=xnew;

        x(s)=x(s)+sum(xnew);
    end

end

r=r/(N*np*dt);       
x=x/(N*np);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract the results for the fit for r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r0=mean(r);
X=r-r0;

X1=(2*dt/period)*sum(X.*cos(2*pi*tp/period));
X2=(2*dt/period)*sum(X.*sin(2*pi*tp/period));

ar=sqrt(X1^2 + X2^2);

% put the angle between -180 and 180
shift=-((X2<0)&(X1>0) - (X2<0)&(X1<0))*pi;
thr=atan2(X1,X2);

rfit=r0 + ar*sin(2*pi*tp/period + thr);

ar=ar./ge1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract the results for the fit for x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0=mean(x);
X=x-x0;

X1=(2*dt/period)*sum(X.*cos(2*pi*tp/period));
X2=(2*dt/period)*sum(X.*sin(2*pi*tp/period));

ax=sqrt(X1^2 + X2^2);

% put the angle between -180 and 180
shift=-((X2<0)&(X1>0) - (X2<0)&(X1<0))*pi;
thx=atan2(X1,X2);

xfit=x0 + ax*sin(2*pi*tp/period + thx);
end