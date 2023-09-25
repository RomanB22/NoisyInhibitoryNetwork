clear all
tau = 10;
a=0.1; b=0.26;
K=1000; % useful for converting rates from kHz to Hz

% parameters for the EIF model
Vre=-60; Vth=30;   % spike currents, threshold and reset
taud=6; taul=1; taur=1;

sig0Sim = [0.51, 1.01:1.0:10.01];

for GoalFreq = [17,30]
    IDX = 1;    
    Esstar=zeros(size(sig0Sim));
    fstar=zeros(size(sig0Sim));
    
    for l = 1:length(sig0Sim)
        File = load(sprintf('FreqResp%dHz/V2Sigma%1.2f.mat',GoalFreq,sig0Sim(l)));
        index=1:size(File.M,1);
        freq=File.M(index,1);
        AE=abs(File.M(index,2));
        phE=medfilt1(File.M(index,3));
               
        w=(freq*2*pi);
        ph2=atan(w*taur)+atan(w*taud)+w*taul-pi;          % NB the +pi is because it is inhibitory !!
        
        [fstarAux,~,k1Aux1,~]=intersectionsCurves(freq,phE,freq,ph2);
        k1Aux = floor(k1Aux1);
        
        EsstarAux = zeros(1,length(k1Aux));
        
        for jj=1:length(k1Aux)
            k1=k1Aux(jj);
            k2=k1+1;
            AEstarAux=AE(k1) + (2*pi*fstarAux(jj)-w(k1))*(AE(k2)-AE(k1))/(w(k2)-w(k1));
            EsstarAux(jj) = -sqrt(1+(2*pi*fstarAux(jj)*taud)^2+(2*pi*fstarAux(jj)*taur)^2)/(tau*AEstarAux);
        end
        
        fstar2(:,1:length(fstarAux))=fstarAux';
        Esstar2(:,1:length(EsstarAux))=EsstarAux;
        
        [Esstar(IDX),Index]=min(abs(Esstar2(:,:)));
        %     [Esstar(IDX),Index]=max(abs(Esstar2(:,:)));
        fstar(IDX)=fstar2(:,Index);
%             Esstar(IDX,1:length(EsstarAux))=EsstarAux;
             fstar(IDX)=fstar2(:,Index);
        IDX=IDX+1;
    end
    M = [sig0Sim;Esstar]'
    figure()
    loglog(sig0Sim,Esstar,'o-')
    ylim([1,153])
    xlim([0.1,10.1])
    %figure()
    %loglog(sig0Sim,fstar,'o-')
end