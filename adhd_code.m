% The code is suitable for our paper "Abnormal hemispheric asymmetry of
% both brain function and structure in attention deficit/hyperactivity
% disorder: a meta-analysis of individual participant data".
% 2021-01-15

all_phe={};  % Demographics of seven sites 
all_phe{1,1}=phekki;
all_phe{1,2}=pheneuro;
all_phe{1,3}=phenyu;
all_phe{1,4}=pheohsu;
all_phe{1,5}=phepek1;
all_phe{1,6}=phepek2;
all_phe{1,7}=phepek3;
site=7;
% age,gender,meanFD ttest2/chi2
for i=1:site
    ad=find(all_phe{1,i}(:,6)==1);
    nn=find(all_phe{1,i}(:,6)==-1);
    [~,P_age(i,1),~,STATS_age{i,1}] = ttest2(all_phe{1,i}(ad,4),all_phe{1,i}(nn,4)); % age ttest2
    group=all_phe{1,i}(:,6);
    gender=all_phe{1,i}(:,3);
    [~,chi2_p(i,1),chi2_p(i,2)] = crosstab(group,gender); % gender chi2test
    [~,P_mfd(i,1),~,STATS_mfd{i,1}] = ttest2(all_phe{1,i}(ad,12),all_phe{1,i}(nn,12)); % meanFD ttest2
end
 %%%%%%%%%% gmv %%%%%%%%%
for i=1:site
    % Asymmetry index calculation
    asy_ad{1,i}=(2*(adgmv{1,i}(:,1:2:111)-adgmv{1,i}(:,2:2:112)))./(adgmv{1,i}(:,1:2:111)+adgmv{1,i}(:,2:2:112));
    asy_nn{1,i}=(2*(nngmv{1,i}(:,1:2:111)-nngmv{1,i}(:,2:2:112)))./(2*(nngmv{1,i}(:,1:2:111)+nngmv{1,i}(:,2:2:112)));
    asy{1,i}=[asy_ad{1,i};asy_nn{1,i}];
    % Within-dataset analysis
    X2=[all_phe{1,i}(:,6),all_phe{1,i}(:,3),all_phe{1,i}(:,4),(all_phe{1,i}(:,4)).^2];  ]
    num=56;
    for j=1:num
        S{j,1} = regstats(asy{1,i}(:,j), X2 ,'linear',{'tstat'}) ;
        gmv7betasetp{1,1}(j,i)=S{j,1}.tstat.beta(2);%beta
        gmv7betasetp{1,2}(j,i)=S{j,1}.tstat.se(2);%se
        gmv7betasetp{1,3}(j,i)=S{j,1}.tstat.t(2);%t
        gmv7betasetp{1,4}(j,i)=S{j,1}.tstat.pval(2);%pval
    end 
end
csvwrite('gmv7p.csv',gmv7betasetp{1,4}) 
%%% P-val was converted to its corresponding Z-score in R. Please see the R Script 'adhd_code_r.R'. 
gmv7ptoz=csvread('gmv7ptoz.csv',1,1); 
a=(gmv7betasetp{1, 3}>0);
gmv7ptoz2=gmv7ptoz;
for i=1:num
    for k=1:site
        if a(i,k)==0
            gmv7ptoz2(i,k)=-gmv7ptoz2(i,k);
        end
    end
end
csvwrite('gmv7ptoz2.csv',gmv7ptoz2) 
% Meta analysis 
%%% The combined Z-score was calculated in R. Please see the R Script 'adhd_code_r.R'.
% Heterogeneity was assessed via the Q statistic
n=[47,66,241,78,90,64,41];
 Q=zeros(56,1);
 zmean=mean((gmv7ptoz2)');
 a=0;
 for k=1:site
     a=a+sqrt(n(k));
 end
 for i=1:num
     for k=1:site
         w(k)=sqrt(n(k))/a;
         Q(i,1)=Q(i,1)+w(k)*((gmv7ptoz2(i,k)-zmean(i))^2);
     end
 end
 csvwrite('Q.csv',Q)
 
% Correlational analysis 
%%% one
masy_gmv_ad=[];
masy_gmv_nn=[];
x1=[];x2=[];
x3=[];x4=[];
for i=1:site
    masy_gmv_ad=[masy_gmv_ad;mean((abs(asy_ad{1,i}))')];
    masy_gmv_nn=[masy_gmv_nn;mean((abs(asy_nn{1,i}))')];
    ad=find(all_phe{1,i}(:,6)==1);
    nn=find(all_phe{1,i}(:,6)==-1);
    x1=[x1;all_phe{1,i}(ad,12)]; % meanFD
    x2=[x2;all_phe{1,i}(ad,13)]; % site
    x3=[x3;all_phe{1,i}(nn,12)]; 
    x4=[x4;all_phe{1,i}(nn,13)];
end 
masy_ad=[masy_gmv_ad,masy_falff_ad,masy_reho_ad];
masy_nn=[masy_gmv_nn,masy_falff_nn,masy_reho_nn];
%%% regressed effects of the site and head motion
for i=1:3
    ad_STATS{1,i} = regstats(masy_ad(:,i),[x1,x2],'linear','r' );
    nn_STATS{1,i} = regstats(masy_nn(:,i),[x3,x4],'linear','r' );
    masy_ad_r(:,i)=ad_STATS{1,i}.r;
    masy_nn_r(:,i)=nn_STATS{1,i}.r;
end
[adr,adp]=corr(masy_ad_r)
[nnr,nnp]=corr(masy_nn_r)
%%% two
phe5=[];
asy_commonROI=[];
x11=[]; x22=[];
for i=1:5
    a=(all_phe{1,i}(:,7:11)-repmat(mean(all_phe{1,i}(:,7:11)),size(all_phe{1,i}(:,7:11),1),1))./repmat(std(all_phe{1,i}(:,7:11)),size(all_phe{1,i}(:,7:11),1),1);
    phe5=[phe5;a]; 
    asy_commonROI=[asy_commonROI;asy{1,i}(:,commonROI)];
    x11=[x11;all_phe{1,i}(:,12)]; 
    x22=[x22;all_phe{1,i}(:,13)];
end
for i=1:length(commonROI)
    STAT{1,i} = regstats(asy_commonROI(:,i),[x11,x22],'linear','r' );
    asy_commonROI_r(:,i)=STAT{1,i}.r;
end
[r,p]=corr(asy_commonROI_r,phe5)

 
 
 