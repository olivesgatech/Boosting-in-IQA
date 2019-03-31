numExp=100;


%%LOAD MAT FILES
% Load Ground Truth
temp=load('dmos_live.mat');
gt_L=temp.y;
temp=load('Multi_GT.mat');
gt_M=temp.y;
temp=load('TID_GT.mat');
gt_T=temp.groundTruth';

temp=load('unique_live.mat');
unique_L=temp.metricRes;
temp=load('unique_multi.mat');
unique_M=temp.temp;
temp=load('unique_tid.mat');
unique_T=temp.temp;
%
cd 'data_L'
names=dir; 
lenN=length(names);
lenGT=length(gt_L);
iqa_L=zeros(lenN,lenGT);
metInd=[5:12,3,4];
for kk=1:10
    temp=load(names(metInd(kk)).name);
    iqa_L(kk,:)=temp.metricRes;
%     metInd(kk)
end
iqa_L(lenN-1,:)=unique_L;
cd ..
%
cd 'data_M'
names=dir; 
lenN=length(names);
lenGT=length(gt_M);
iqa_M=zeros(lenN,lenGT);
metInd=[5:12,3,4];
for kk=1:10
    temp=load(names(metInd(kk)).name);
    iqa_M(kk,:)=temp.temp;
%     metInd(kk)
end
iqa_M(lenN-1,:)=unique_M;
cd ..

cd 'data_T'
names=dir; 
lenN=length(names);
lenGT=length(gt_T);
iqa_T=zeros(lenN,lenGT);
metInd=[5:12,3,4];
for kk=1:10
    temp=load(names(metInd(kk)).name);
    iqa_T(kk,:)=temp.temp;
%     metInd(kk)
end
iqa_T(lenN-1,:)=unique_T;
cd ..
%% REGRESS QUALITY ESTIAMTES
iqa_L=regressMethods(iqa_L,gt_L);
iqa_M=regressMethods(iqa_M,gt_M);
iqa_T=regressMethods(iqa_T,gt_T);

%% PART 1
nMethod=12;
corr_L_SVR=zeros(nMethod,lenN,numExp);
corr_M_SVR=zeros(nMethod,lenN,numExp);
corr_T_SVR=zeros(nMethod,lenN,numExp);
P_corr_L_SVR=zeros(nMethod,lenN,numExp);
P_corr_M_SVR=zeros(nMethod,lenN,numExp);
P_corr_T_SVR=zeros(nMethod,lenN,numExp);
mse_L_SVR=zeros(nMethod,lenN,numExp);
mse_M_SVR=zeros(nMethod,lenN,numExp);
mse_T_SVR=zeros(nMethod,lenN,numExp);

K=5;
%%SVR
for pp=1:nMethod
    if (pp<12)
    f_L=iqa_L(pp,:)';
    f_M=iqa_M(pp,:)'; 
    f_T=iqa_T(pp,:)';    
    
    else
    f_L=iqa_L(1:end-1,:)';
    f_M=iqa_M(1:end-1,:)';
    f_T=iqa_T(1:end-1,:)';
    
    end
 
%     f_L_t=f_L';
%     gt_L_t=gt_L';

for mm=1:numExp
    MethodNo=pp
    ExperimentNo=mm
    % LIVE
    len=length(gt_L);
    valInd=crossvalind('Kfold', len, K);
    counter_TR=1;
    counter_TE=1;
    test_f_L=[];
    train_f_L=[];
    train_gt_L=[];
    test_gt_L=[];
   if pp==12      
         for ii=1:len
            if(valInd(ii)==5)
                test_f_L(:,counter_TE)=f_L(ii,:);
                test_gt_L(counter_TE)=gt_L(ii);
                counter_TE=counter_TE+1;
            else
                train_f_L(:,counter_TR)=f_L(ii,:);
                train_gt_L(counter_TR)=gt_L(ii);
                counter_TR=counter_TR+1;            
            end
        end
       
   else       
        for ii=1:len
            if(valInd(ii)==5)
                test_f_L(counter_TE)=f_L(ii);
                test_gt_L(counter_TE)=gt_L(ii);
                counter_TE=counter_TE+1;
            else
                train_f_L(counter_TR)=f_L(ii);
                train_gt_L(counter_TR)=gt_L(ii);
                counter_TR=counter_TR+1;            
            end
        end
   end
    
    testInd=find(valInd==5);
    
    mdl_L=fitrsvm(train_f_L',train_gt_L);
    est_L=predict(mdl_L,test_f_L');

    test_iqa_L=iqa_L(:,testInd);
    test_iqa_L(end,:)=est_L';


for nn=1:lenN
   corr_L_SVR(pp,nn,mm)=abs(corr(test_iqa_L(nn,:)',test_gt_L','Type','Spearman'));
   mse_L_SVR(pp,nn,mm)=abs(mse_1D(test_iqa_L(nn,:)',test_gt_L'));
   P_corr_L_SVR(pp,nn,mm)=abs(corr(test_iqa_L(nn,:)',test_gt_L','Type','Pearson'));

end

%M
    len=length(gt_M);
    valInd=crossvalind('Kfold', len, K);
    counter_TR=1;
    counter_TE=1;
    test_f_M=[];
    train_f_M=[];
    train_gt_M=[];
    test_gt_M=[];
   if pp==12      
         for ii=1:len
            if(valInd(ii)==5)
                test_f_M(:,counter_TE)=f_M(ii,:);
                test_gt_M(counter_TE)=gt_M(ii);
                counter_TE=counter_TE+1;
            else
                train_f_M(:,counter_TR)=f_M(ii,:);
                train_gt_M(counter_TR)=gt_M(ii);
                counter_TR=counter_TR+1;            
            end
        end
       
   else       
        for ii=1:len
            if(valInd(ii)==5)
                test_f_M(counter_TE)=f_M(ii);
                test_gt_M(counter_TE)=gt_M(ii);
                counter_TE=counter_TE+1;
            else
                train_f_M(counter_TR)=f_M(ii);
                train_gt_M(counter_TR)=gt_M(ii);
                counter_TR=counter_TR+1;            
            end
        end
   end
    
    testInd=find(valInd==5);
    
    mdl_M=fitrsvm(train_f_M',train_gt_M');
    est_M=predict(mdl_M,test_f_M');

test_iqa_M=iqa_M(:,testInd);
test_iqa_M(end,:)=est_M';

for nn=1:lenN
   corr_M_SVR(pp,nn,mm)=abs(corr(test_iqa_M(nn,:)',test_gt_M','Type','Spearman'));
   mse_M_SVR(pp,nn,mm)=abs(mse_1D(test_iqa_M(nn,:)',test_gt_M'));   
   P_corr_M_SVR(pp,nn,mm)=abs(corr(test_iqa_M(nn,:)',test_gt_M','Type','Pearson'));
   
end

%T
    len=length(gt_T);
    valInd=crossvalind('Kfold', len, K);
    counter_TR=1;
    counter_TE=1;
    test_f_T=[];
    train_f_T=[];
    train_gt_T=[];
    test_gt_T=[];
   if pp==12      
         for ii=1:len
            if(valInd(ii)==5)
                test_f_T(:,counter_TE)=f_T(ii,:);
                test_gt_T(counter_TE)=gt_T(ii);
                counter_TE=counter_TE+1;
            else
                train_f_T(:,counter_TR)=f_T(ii,:);
                train_gt_T(counter_TR)=gt_T(ii);
                counter_TR=counter_TR+1;            
            end
        end
       
   else       
        for ii=1:len
            if(valInd(ii)==5)
                test_f_T(counter_TE)=f_T(ii);
                test_gt_T(counter_TE)=gt_T(ii);
                counter_TE=counter_TE+1;
            else
                train_f_T(counter_TR)=f_T(ii);
                train_gt_T(counter_TR)=gt_T(ii);
                counter_TR=counter_TR+1;            
            end
        end
   end
    
    testInd=find(valInd==5);
    
    mdl_T=fitrsvm(train_f_T',train_gt_T');
    est_T=predict(mdl_T,test_f_T');

test_iqa_T=iqa_T(:,testInd);
test_iqa_T(end,:)=est_T';

for nn=1:lenN
   corr_T_SVR(pp,nn,mm)=abs(corr(test_iqa_T(nn,:)',test_gt_T','Type','Spearman'));
   mse_T_SVR(pp,nn,mm)=abs(mse_1D(test_iqa_T(nn,:)',test_gt_T'));   
   P_corr_T_SVR(pp,nn,mm)=abs(corr(test_iqa_T(nn,:)',test_gt_T','Type','Pearson'));   
end
end
end



% %not trained - MSE
% squeeze(mean(squeeze(mean(mse_L_SVR(:,1:(end-1),:),1)),2))'
% squeeze(mean(squeeze(mean(mse_M_SVR(:,1:(end-1),:),1)),2))'
% squeeze(mean(squeeze(mean(mse_T_SVR(:,1:(end-1),:),1)),2))'
% %not trained - P CC
% squeeze(mean(squeeze(mean(P_corr_L_SVR(:,1:(end-1),:),1)),2))'
% squeeze(mean(squeeze(mean(P_corr_M_SVR(:,1:(end-1),:),1)),2))'
% squeeze(mean(squeeze(mean(P_corr_T_SVR(:,1:(end-1),:),1)),2))'
% %not trained - CC
% squeeze(mean(squeeze(mean(corr_L_SVR(:,1:(end-1),:),1)),2))'
% squeeze(mean(squeeze(mean(corr_M_SVR(:,1:(end-1),:),1)),2))'
% squeeze(mean(squeeze(mean(corr_T_SVR(:,1:(end-1),:),1)),2))'
% 

% trained - MSE
squeeze(mean(mse_L_SVR(:,end,:),3))'
squeeze(mean(mse_M_SVR(:,end,:),3))'
squeeze(mean(mse_T_SVR(:,end,:),3))'
% trained - P CC
squeeze(mean(P_corr_L_SVR(:,end,:),3))'
squeeze(mean(P_corr_M_SVR(:,end,:),3))'
squeeze(mean(P_corr_T_SVR(:,end,:),3))'
% trained - CC
squeeze(mean(corr_L_SVR(:,end,:),3))'
squeeze(mean(corr_M_SVR(:,end,:),3))'
squeeze(mean(corr_T_SVR(:,end,:),3))'
% %% Part 2
% 
% load('orders.mat');
% % o1=fliplr(o1);
% % o2=fliplr(o2);
% % o3=fliplr(o3);
% % o4=fliplr(o4);
% % o5=fliplr(o5);
% % o6=fliplr(o6);
% % o7=fliplr(o7);
% % o8=fliplr(o8);
% % o9=fliplr(o9);
% 
% ord=vertcat(o1,o2,o3,o4,o5,o6,o7,o8,o9);
% 
% nMethod=3;
% lenN=11;
% K=5;
% corr_L_SVR=zeros(nMethod,lenN,numExp);
% corr_M_SVR=zeros(nMethod,lenN,numExp);
% corr_T_SVR=zeros(nMethod,lenN,numExp);
% 
% P_corr_L_SVR=zeros(nMethod,lenN,numExp);
% P_corr_M_SVR=zeros(nMethod,lenN,numExp);
% P_corr_T_SVR=zeros(nMethod,lenN,numExp);
% 
% mse_L_SVR=zeros(nMethod,lenN,numExp);
% mse_M_SVR=zeros(nMethod,lenN,numExp);
% mse_T_SVR=zeros(nMethod,lenN,numExp);
% 
% for pp=1:nMethod
% f_L=[];
% f_M=[];
% f_T=[];
% 
%     for kk=1:lenN
%     ord1=(pp-1)*3+1;
%     ord2=(pp-1)*3+2;
%     ord3=pp*3;
%     f_L=horzcat(f_L,iqa_L(ord(ord1,kk),:)');  
%     f_M=horzcat(f_M,iqa_M(ord(ord2,kk),:)');  
%     f_T=horzcat(f_T,iqa_T(ord(ord3,kk),:)');  
% 
%         for mm=1:numExp
%         
%             ValMetricNo=pp
%             ExperimentNo=mm
%             
%             %LIVE
%             f_L_t=f_L';
%             gt_L_t=gt_L';
% 
%             len=length(gt_L);
%             valInd=crossvalind('Kfold', len, K);
%             counter_TR=1;
%             counter_TE=1;
%             test_f_L=[];
%             train_f_L=[];
%             train_gt_L=[];
%             test_gt_L=[];
%             
%             for ii=1:len
%                 if(valInd(ii)==5)
%                     test_f_L(:,counter_TE)=f_L(ii,:);
%                     test_gt_L(counter_TE)=gt_L(ii);
%                     counter_TE=counter_TE+1;
%                 else
%                     train_f_L(:,counter_TR)=f_L(ii,:);
%                     train_gt_L(counter_TR)=gt_L(ii);
%                     counter_TR=counter_TR+1;            
%                 end
%             end
%             
%             testInd=find(valInd==5);
%               
%             mdl_L=fitrsvm(train_f_L',train_gt_L);
%             est_L=predict(mdl_L,test_f_L');
%            
%    corr_L_SVR(pp,kk,mm)=abs(corr(est_L,test_gt_L','Type','Spearman'));
%    mse_L_SVR(pp,kk,mm)=abs(mse_1D(est_L,test_gt_L'));
%    P_corr_L_SVR(pp,kk,mm)=abs(corr(est_L,test_gt_L','Type','Pearson'));
%             %% MULTI
%             f_M_t=f_M';
%             gt_M_t=gt_M';
% 
%             len=length(gt_M);
%             valInd=crossvalind('Kfold', len, K);
%             counter_TR=1;
%             counter_TE=1;
%             test_f_M=[];
%             train_f_M=[];
%             train_gt_M=[];
%             test_gt_M=[];
%             
%             for ii=1:len
% 
%                 if(valInd(ii)==5)
%                     test_f_M(:,counter_TE)=f_M(ii,:);
%                     test_gt_M(counter_TE)=gt_M(ii);
%                     counter_TE=counter_TE+1;
%                 else
%                     train_f_M(:,counter_TR)=f_M(ii,:);
%                     train_gt_M(counter_TR)=gt_M(ii);
%                  counter_TR=counter_TR+1;            
%                 end
%             end
%             
%             testInd=find(valInd==5);
%               
%             mdl_M=fitrsvm(train_f_M',train_gt_M);
%             est_M=predict(mdl_M,test_f_M');
%            
%    corr_M_SVR(pp,kk,mm)=abs(corr(est_M,test_gt_M','Type','Spearman'));
%    mse_M_SVR(pp,kk,mm)=abs(mse_1D(est_M,test_gt_M'));
%    P_corr_M_SVR(pp,kk,mm)=abs(corr(est_M,test_gt_M','Type','Pearson'));
%             
%            % TID
%             f_T_t=f_T';
%             gt_T_t=gt_T';
% 
%             len=length(gt_T);
%             valInd=crossvalind('Kfold', len, K);
%             counter_TR=1;
%             counter_TE=1;
%             test_f_T=[];
%             train_f_T=[];
%             train_gt_T=[];
%             test_gt_T=[];
%             
%             for ii=1:len
% 
%                 if(valInd(ii)==5)
%                     test_f_T(:,counter_TE)=f_T(ii,:);
%                     test_gt_T(counter_TE)=gt_T(ii);
%                     counter_TE=counter_TE+1;
%                 else
%                     train_f_T(:,counter_TR)=f_T(ii,:);
%                     train_gt_T(counter_TR)=gt_T(ii);
%                  counter_TR=counter_TR+1;            
%                 end
%             end
%             
%             testInd=find(valInd==5);
%               
%             mdl_T=fitrsvm(train_f_T',train_gt_T);
%             est_T=predict(mdl_T,test_f_T');
%            
%    corr_T_SVR(pp,kk,mm)=abs(corr(est_T,test_gt_T','Type','Spearman'));
%    mse_T_SVR(pp,kk,mm)=abs(mse_1D(est_T,test_gt_T'));
%    P_corr_T_SVR(pp,kk,mm)=abs(corr(est_T,test_gt_T','Type','Pearson'));  
%             
%             %%
%             
%         end
%     end
% end
% 
% 
% 
% 
% 
% % trained - MSE
% mean_mse_L=squeeze(mean(mse_L_SVR(1,:,:),3))'
% mean_mse_M=squeeze(mean(mse_M_SVR(1,:,:),3))'
% mean_mse_T=squeeze(mean(mse_T_SVR(1,:,:),3))'
% 
% % trained -P CC
% mean_P_L=squeeze(mean(P_corr_L_SVR(2,:,:),3))'
% mean_P_M=squeeze(mean(P_corr_M_SVR(2,:,:),3))'
% mean_P_T=squeeze(mean(P_corr_T_SVR(2,:,:),3))'
% 
% % trained - CC
% mean_S_L=squeeze(mean(corr_L_SVR(3,:,:),3))'
% mean_S_M=squeeze(mean(corr_M_SVR(3,:,:),3))'
% mean_S_T=squeeze(mean(corr_T_SVR(3,:,:),3))'
% 
% % trained - MSE
% std_mse_L=squeeze(std(mse_L_SVR(1,:,:),0,3))'
% std_mse_M=squeeze(std(mse_M_SVR(1,:,:),0,3))'
% std_mse_T=squeeze(std(mse_T_SVR(1,:,:),0,3))'
% 
% % trained -P CC
% std_P_L=squeeze(std(P_corr_L_SVR(2,:,:),0,3))'
% std_P_M=squeeze(std(P_corr_M_SVR(2,:,:),0,3))'
% std_P_T=squeeze(std(P_corr_T_SVR(2,:,:),0,3))'
% 
% % trained - CC
% std_S_L=squeeze(std(corr_L_SVR(3,:,:),0,3))'
% std_S_M=squeeze(std(corr_M_SVR(3,:,:),0,3))'
% std_S_T=squeeze(std(corr_T_SVR(3,:,:),0,3))'
% 
% 
% 
% 
% 
% 
% %% Figures
% scType='png';
% len=length(mean_mse_L);
% 
% %MSE - L
% figure
% hold on
% bar(1:len,mean_mse_L);
% ax=gca;
% ax.FontSize=20;
% errorbar(1:len,mean_mse_L,std_mse_L,'.','Linewidth',4,'color','r')
% xlabel('Number of fused methods','fontsize',30,'fontweight','bold');
% ylabel('RMSE','fontsize',30,'fontweight','bold');
% axis([0,12,5,13])
% saveas(gcf,'SVR_mse_L',scType)
% 
% %MSE - M
% figure
% hold on
% bar(1:len,mean_mse_M);
% ax=gca;
% ax.FontSize=20;
% errorbar(1:len,mean_mse_M,std_mse_M,'.','Linewidth',4,'color','r')
% xlabel('Number of fused methods','fontsize',30,'fontweight','bold'); 
% ylabel('RMSE','fontsize',30,'fontweight','bold');
% axis([0,12,6,21])
% saveas(gcf,'SVR_mse_M',scType)
% 
% %MSE - T
% figure
% hold on
% bar(1:len,mean_mse_T);
% ax=gca;
% ax.FontSize=20;
% errorbar(1:len,mean_mse_T,std_mse_T,'.','Linewidth',4,'color','r')
% xlabel('Number of fused methods','fontsize',30,'fontweight','bold'); 
% ylabel('RMSE','fontsize',30,'fontweight','bold');
% axis([0,12,0,1.4])
% saveas(gcf,'SVR_mse_T',scType)
% 
% 
% R1=mean_P_L(1);
% alpha=0.05;
% nSamp=length(gt_L);
% sigmaZ=sqrt(1/(nSamp-3));
% z1=0.5*log((1+R1)/(1-R1));
% for ii=R1:0.01:1.00
% z2=0.5*log((1+ii)/(1-ii));
% zn=abs(z1-z2)/sqrt(2*sigmaZ^2);
%     if (zn > tinv(1-alpha/2,nSamp-1))
%     break
%     else
%     end
% end
% %P - L
% figure
% hold on
% bar(1:len,mean_P_L);
% ax=gca;
% ax.FontSize=20;
% errorbar(1:len,mean_P_L,std_P_L,'.','Linewidth',4,'color','r')
% hold on
% hline(ii,'k-x','SS')
% xlabel('Number of fused methods','fontsize',30,'fontweight','bold');
% ylabel('PLCC','fontsize',30,'fontweight','bold');
% axis([0,12,0.85,0.98])
% saveas(gcf,'SVR_P_L',scType)
% 
% 
% R1=mean_P_M(1);
% alpha=0.05;
% nSamp=length(gt_M);
% sigmaZ=sqrt(1/(nSamp-3));
% z1=0.5*log((1+R1)/(1-R1));
% for ii=R1:0.01:1.00
% z2=0.5*log((1+ii)/(1-ii));
% zn=abs(z1-z2)/sqrt(2*sigmaZ^2);
%     if (zn > tinv(1-alpha/2,nSamp-1))
%     break
%     else
%     end
% end
% %P - M
% figure
% hold on
% bar(1:len,mean_P_M);
% ax=gca;
% ax.FontSize=20;
% errorbar(1:len,mean_P_M,std_P_M,'.','Linewidth',4,'color','r')
% hold on
% hline(ii,'k-x','SS')
% xlabel('Number of fused methods','fontsize',30,'fontweight','bold'); 
% ylabel('PLCC','fontsize',30,'fontweight','bold');
% axis([0,12,0.25,1.00])
% saveas(gcf,'SVR_P_M',scType)
% 
% R1=mean_P_T(1);
% alpha=0.05;
% nSamp=length(gt_T);
% sigmaZ=sqrt(1/(nSamp-3));
% z1=0.5*log((1+R1)/(1-R1));
% for ii=R1:0.01:1.00
% z2=0.5*log((1+ii)/(1-ii));
% zn=abs(z1-z2)/sqrt(2*sigmaZ^2);
%     if (zn > tinv(1-alpha/2,nSamp-1))
%     break
%     else
%     end
% end
% %P - M
% figure
% hold on
% bar(1:len,mean_P_T);
% ax=gca;
% ax.FontSize=20;
% errorbar(1:len,mean_P_T,std_P_T,'.','Linewidth',4,'color','r')
% hold on
% hline(ii,'k-x','SS')
% xlabel('Number of fused methods','fontsize',30,'fontweight','bold'); 
% ylabel('PLCC','fontsize',30,'fontweight','bold');
% axis([0,12,0.1,1.00])
% saveas(gcf,'SVR_P_T',scType)
% 
% 
% 
% 
% 
% R1=mean_S_L(1);
% alpha=0.05;
% nSamp=length(gt_L);
% sigmaZ=sqrt(1/(nSamp-3));
% z1=0.5*log((1+R1)/(1-R1));
% for ii=R1:0.01:1.00
% z2=0.5*log((1+ii)/(1-ii));
% zn=abs(z1-z2)/sqrt(2*sigmaZ^2);
%     if (zn > tinv(1-alpha/2,nSamp-1))
%     break
%     else
%     end
% end
% %S - L
% figure
% hold on
% bar(1:len,mean_S_L);
% ax=gca;
% ax.FontSize=20;
% errorbar(1:len,mean_S_L,std_S_L,'.','Linewidth',4,'color','r')
% hold on
% hline(ii,'k-x','SS')
% xlabel('Number of fused methods','fontsize',30,'fontweight','bold');
% ylabel('SRCC','fontsize',30,'fontweight','bold');
% axis([0,12,0.88,0.965])
% saveas(gcf,'SVR_S_L',scType)
% 
% 
% R1=mean_S_M(1);
% alpha=0.05;
% nSamp=length(gt_M);
% sigmaZ=sqrt(1/(nSamp-3));
% z1=0.5*log((1+R1)/(1-R1));
% for ii=R1:0.01:1.00
% z2=0.5*log((1+ii)/(1-ii));
% zn=abs(z1-z2)/sqrt(2*sigmaZ^2);
%     if (zn > tinv(1-alpha/2,nSamp-1))
%     break
%     else
%     end
% end
% %S - M
% figure
% hold on
% bar(1:len,mean_S_M);
% ax=gca;
% ax.FontSize=20;
% errorbar(1:len,mean_S_M,std_S_M,'.','Linewidth',4,'color','r')
% hold on
% hline(ii,'k-x','SS')
% xlabel('Number of fused methods','fontsize',30,'fontweight','bold'); 
% ylabel('SRCC','fontsize',30,'fontweight','bold');
% axis([0,12,0.52,1.00])
% saveas(gcf,'SVR_S_M',scType)
% 
% 
% 
% 
% R1=mean_S_T(1);
% alpha=0.05;
% nSamp=length(gt_T);
% sigmaZ=sqrt(1/(nSamp-3));
% z1=0.5*log((1+R1)/(1-R1));
% for ii=R1:0.01:1.00
% z2=0.5*log((1+ii)/(1-ii));
% zn=abs(z1-z2)/sqrt(2*sigmaZ^2);
%     if (zn > tinv(1-alpha/2,nSamp-1))
%     break
%     else
%     end
% end
% %S - M
% figure
% hold on
% bar(1:len,mean_S_T);
% ax=gca;
% ax.FontSize=20;
% errorbar(1:len,mean_S_T,std_S_T,'.','Linewidth',4,'color','r')
% hold on
% hline(ii,'k-x','SS')
% xlabel('Number of fused methods','fontsize',30,'fontweight','bold'); 
% ylabel('SRCC','fontsize',30,'fontweight','bold');
% axis([0,12,0.52,0.95])
% saveas(gcf,'SVR_S_T',scType)
% 














