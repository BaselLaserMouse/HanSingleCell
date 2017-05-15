%% consider the set of neurons that project to two areas only. in this set are there certain classes/motifs
% that are under or overrepresented relative to what we would expect if the
% neurons in this set were distributed according to the product of the
% first order innervation probabilities?
function [motif,pval_min,significance,overorunder,observed,expected,p_all,neuronnum]=higherorder(data,cutoff,alpha)



%binarize input matrix
data_binary_tmp=double(data>=cutoff);
data_binary=data_binary_tmp(sum(data_binary_tmp,2)~=0,:);

targetnum=size(data_binary,2);

%uniq
[data_binary_uniq,ia,ib]=unique(data_binary,'rows');
data_counts=accumarray(ib,1);



%get first order probabilities


neuronnum=findntotal(data,cutoff);

p=sum(data_binary,1)./neuronnum;



%construct motif matrix and associated expected probabilities from first
%order stats

%bifrucations
motif_tmp=zeros(targetnum^2,targetnum);
p_expected_tmp=zeros(targetnum^2,1);
counter=1;
for target1=1:targetnum;
    for target2=1:targetnum;
        motif_tmp(counter,target1)=1;
        motif_tmp(counter,target2)=1;
        p_expected_tmp(counter)=prod(p(logical(motif_tmp(counter,:))))*prod(1-p(~logical(motif_tmp(counter,:))));
        counter=counter+1;
    end
end
%remove doublecounts
[motif_bi,ia,ib]=unique(motif_tmp(sum(motif_tmp,2)==2,:),'rows');
p_expected_tmp2=p_expected_tmp(sum(motif_tmp,2)==2);
p_expected_tmp3_bi=p_expected_tmp2(ia);


%trifrucations
motif_tmp=zeros(targetnum^3,targetnum);
p_expected_tmp=zeros(targetnum^3,1);
counter=1;
for target1=1:targetnum;
    for target2=1:targetnum;
        for target3=1:targetnum;
            motif_tmp(counter,target1)=1;
            motif_tmp(counter,target2)=1;
            motif_tmp(counter,target3)=1;
            p_expected_tmp(counter)=prod(p(logical(motif_tmp(counter,:))))*prod(1-p(~logical(motif_tmp(counter,:))));            
            counter=counter+1;
        end
    end
end
%remove doublecounts
[motif_tri,ia,ib]=unique(motif_tmp(sum(motif_tmp,2)==3,:),'rows');
p_expected_tmp2=p_expected_tmp(sum(motif_tmp,2)==3);
p_expected_tmp3_tri=p_expected_tmp2(ia);

% quadfrucations
motif_tmp=zeros(targetnum^4,targetnum);
p_expected_tmp=zeros(targetnum^4,1);
counter=1;
for target1=1:targetnum;
    for target2=1:targetnum;
        for target3=1:targetnum;
            for target4=1:targetnum;
                motif_tmp(counter,target1)=1;
                motif_tmp(counter,target2)=1;
                motif_tmp(counter,target3)=1;
                motif_tmp(counter,target4)=1;
                p_expected_tmp(counter)=prod(p(logical(motif_tmp(counter,:))))*prod(1-p(~logical(motif_tmp(counter,:))));                
                counter=counter+1;
            end
        end
    end
end
%remove doublecounts
[motif_quad,ia,ib]=unique(motif_tmp(sum(motif_tmp,2)==4,:),'rows');
p_expected_tmp2=p_expected_tmp(sum(motif_tmp,2)==4);
p_expected_tmp3_quad=p_expected_tmp2(ia);

%combine all motifs and probabilities
motif=[motif_bi;motif_tri;motif_quad];
p_all=[p_expected_tmp3_bi;p_expected_tmp3_tri;p_expected_tmp3_quad];



%count up the number of observed 
observed=zeros(size(p_all));

[d,loc]=ismember(motif,data_binary_uniq,'rows');
observed(d)=data_counts(loc(loc~=0));

%calculate expected
expected=neuronnum.*p_all;

%calculate pvalues


% pvalue=binocdf(observed,bifrucnum,p_expected);
pvalue=binocdf(observed,ceil(neuronnum),p_all);

%decide on significance using bonferoni correction
motifnum=size(motif,1);
alpha_adj=alpha/motifnum;
overorunder=zeros(size(p_all));
overorunder(pvalue<(alpha_adj/2))=-1;
overorunder(pvalue>(1-alpha_adj/2))=1;


significance=overorunder~=0;

% make plots

%plot observed vs expected
figure;bar([observed,expected]);
xlabel('motif')
ylabel('count')
title('observed and expected counts of polyfrucations')
legend('observed','expected');

%plot volcano plot
pval_min=nan*ones(size(pvalue));
pval_min(pvalue<0.5)=pvalue(pvalue<0.5);
pval_min(pvalue>0.5)=1-pvalue(pvalue>0.5);
figure;mavolcanoplot(observed,expected,pval_min,'LogTrans',1,'PlotOnly',1,'PCutoff',alpha_adj/2)


        
