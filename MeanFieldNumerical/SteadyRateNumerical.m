function r0sim = SteadyRateNumerical(Es0,sig0Sim)

tau = 10;
a=0.1; b=0.26;
K=1000; % useful for converting rates from kHz to Hz

% parameters for the EIF model
EL=-60; Vre=-60; Vth=30;   % spike currents, threshold and reset      

nEs0=length(Es0);

r0sim=zeros(nEs0,1); 
x0sim=zeros(nEs0,1);

for k=1:nEs0
    
    ge0=Es0(k);    
    
    [r0sim(k),x0sim(k)]=simUncNet(ge0,sig0Sim,Vre,Vth,EL,a,b,tau);
    
end