files = dir('FICurve/*.txt');
GoalFreq = [17,30];
sig0Sim = 0.01:0.5:10.01;
sig0Sim = sig0Sim([1:4,21,5:20]);

for i = 1:length(files)
    Matrix = load(strcat(files(i).folder,'/',files(i).name));
    I0 = Matrix(:,1); r0 = Matrix(:,2);
    for GoalFreqAux=GoalFreq
        [I0point,r0point,iout,jout] = intersectionsCurves(I0,r0,I0,GoalFreqAux*ones(size(I0)));
        I0point=min(I0point);
        save(sprintf('I0%dHz/I0Sigma%1.2f.txt',GoalFreqAux,sig0Sim(i)),'I0point','-ascii')
    end
end


