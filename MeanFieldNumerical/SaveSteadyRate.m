Es0=[-30:0.1:10]';
K=1000;

for sig0Sim = 0.01:0.5:10.01
    r0sim = SteadyRateNumerical(Es0,sig0Sim);
    M = [Es0,K*r0sim];
    save(sprintf('FICurve/Sigma%1.2f.txt',sig0Sim),'M','-ascii')
end