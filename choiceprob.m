function [Choiceprobability,hitrate,farate] = choiceprob(X_sig,X_noise,n,o,p,q)
%Function to obtain choice probability value of two distributions. Choice probability (CP) provides information about how well an ideal observer can predict 
%the choice an animal makes from a neuron's discharge rate distribution.
%Typically CPs are meaningful when a stimulus in the external world is
%ambiguous/noisy and a neuron's discharge rate in response to the stimulus is different based on the the animal's choice/judgement about the stimulus.
%Some refs:  Bradley et al. 1987; Britten et al. 1992; Newsome et al. 1989; Parker and Hawken 1985; Skottun et al. 1987; Tolhurst et al. 1983; Vogels and Orban 1991

%  Input: X_sig and X_noise are the list of values of distribution one and
%  two. n,o,p,q are related figure number and subplotting of
%  Receiver-operator characteristic curves (ROC).

%  Ouput: Choiceprobability is CP ranging from 0 to 1. hitrate: percent correct choices indicating the presence of a stimulus feature 
%farate (aka false alarm rate): percent wrong choices indicating the presence of a stimulus feature 
lowfr=min([X_sig,X_noise]);
highfr=max([X_sig,X_noise]);


figure(1)
clf(1)
h_sig = histogram(X_sig,'BinMethod','fd');
min_sig=h_sig.BinWidth;
hold on
h_noise = histogram(X_noise,'BinMethod','fd');
min_noise=h_noise.BinWidth; 

%BW=min(min_sig, min_noise); %keeping same binwidth for sig and noise

BW=0.25; %binwidth

figure(100+n)
%clf(100+n)
subplot(o,p,q)
hold on
h_sig = histogram(X_sig,'BinWidth',BW,'FaceColor','g');
h_noise = histogram(X_noise,'BinWidth',BW,'FaceColor','b');
hold on
legend(sprintf("choice1 n=%d",length(X_sig)),sprintf("choice2 n=%d",length(X_noise)))
xlabel("Response")
ylabel("#Trials")

step=BW;
dcriterion=lowfr:step:highfr;
%area = sum(h.Values)*h.BinWidth;
%histogram(...,'BinLimits',[BMIN,BMAX])
%histogram(...,'BinWidth',BW)
for i=1:length(dcriterion)
    figure(200)
    clf(200)
    if dcriterion(i)<max(X_sig)
       h_hitrate=histogram(X_sig,'BinLimits',[dcriterion(i),max(X_sig)],'BinWidth',BW,'FaceColor','r');
       hold on
       hitrate(i)=sum(h_hitrate.Values)/sum(h_sig.Values);
    else
        hitrate(i)=0;
    end
    if dcriterion(i)<max(X_noise)
       h_farate=histogram(X_noise,'BinLimits',[dcriterion(i),max(X_noise)],'BinWidth',BW,'FaceColor','k');
       farate(i)=sum(h_farate.Values)/sum(h_noise.Values);
    else 
       farate(i)=0;
    end
end

farate=[0 sort(farate)];
hitrate=[0 sort(hitrate)];


%Z = trapz(X,Y)
Choiceprobability =trapz(farate,hitrate);
%dprime = inputArg2;

