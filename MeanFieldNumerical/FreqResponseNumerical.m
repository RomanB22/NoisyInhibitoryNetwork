GoalFreq = [30];

tau = 10;
a=0.1; b=0.26;
K=1000; % useful for converting rates from kHz to Hz
% parameters for the EIF model
EL=-60; Vre=-60; Vth=30;   % spike currents, threshold and reset

lfsim=1/1000; hfsim=250/1000; nfsim=125;
fsim=linspace(lfsim,hfsim,nfsim)';%(exp(linspace(log(lfsim),log(hfsim),nfsim))');

T=50000;
dt1=0.5;
N=20000;

for GoalFreqAux=GoalFreq
    for sig0Sim = 0.01:0.5:10.01
        % general synaptic parameters
        ge0=load(sprintf('I0%dHz/I0Sigma%1.2f.txt',GoalFreqAux,sig0Sim));
        ge1=abs(0.1*ge0);
        
        r0sim=zeros(nfsim,1);
        arsim=zeros(nfsim,1);
        thrsim=zeros(nfsim,1);        
        axsim=zeros(nfsim,1);
        thxsim=zeros(nfsim,1);        
        for k=1:nfsim 
            [r0sim(k),arsim(k),thrsim(k),axsim(k),thxsim(k)]=ge_mod_sim(fsim(k),T,dt1,N,EL,ge0,ge1,Vre,Vth,sig0Sim,a,b,tau);
        end
        M = [fsim,arsim,thrsim,axsim,thxsim];
        save(sprintf('FreqResp%dHz/V2Sigma%1.2f.mat',GoalFreqAux,sig0Sim),'M');
    end    
end
