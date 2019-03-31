% load('table_nn.mat');
% load('table_svr.mat');

%not trained - MSE
S_r1=squeeze(mean(squeeze(mean(mse_L_SVR(:,1:(end-1),:),1)),2))'
S_r2=squeeze(mean(squeeze(mean(mse_M_SVR(:,1:(end-1),:),1)),2))'
S_r3=squeeze(mean(squeeze(mean(mse_T_SVR(:,1:(end-1),:),1)),2))'

%not trained - P CC
S_r4=squeeze(mean(squeeze(mean(P_corr_L_SVR(:,1:(end-1),:),1)),2))'
S_r5=squeeze(mean(squeeze(mean(P_corr_M_SVR(:,1:(end-1),:),1)),2))'
S_r6=squeeze(mean(squeeze(mean(P_corr_T_SVR(:,1:(end-1),:),1)),2))'

%not trained - CC
S_r7=squeeze(mean(squeeze(mean(corr_L_SVR(:,1:(end-1),:),1)),2))'
S_r8=squeeze(mean(squeeze(mean(corr_M_SVR(:,1:(end-1),:),1)),2))'
S_r9=squeeze(mean(squeeze(mean(corr_T_SVR(:,1:(end-1),:),1)),2))'



%not trained - MSE
r1=squeeze(mean(squeeze(mean(mse_L(:,1:(end-1),:),1)),2))'
r2=squeeze(mean(squeeze(mean(mse_M(:,1:(end-1),:),1)),2))'
r3=squeeze(mean(squeeze(mean(mse_T(:,1:(end-1),:),1)),2))'

%not trained - P CC
r4=squeeze(mean(squeeze(mean(P_corr_L(:,1:(end-1),:),1)),2))'
r5=squeeze(mean(squeeze(mean(P_corr_M(:,1:(end-1),:),1)),2))'
r6=squeeze(mean(squeeze(mean(P_corr_T(:,1:(end-1),:),1)),2))'
%not trained - CC
r7=squeeze(mean(squeeze(mean(corr_L(:,1:(end-1),:),1)),2))'
r8=squeeze(mean(squeeze(mean(corr_M(:,1:(end-1),:),1)),2))'
r9=squeeze(mean(squeeze(mean(corr_T(:,1:(end-1),:),1)),2))'


row1=(r1+S_r1)/2
row2=(r2+S_r2)/2
row3=(r3+S_r3)/2
row4=(r4+S_r4)/2
row5=(r5+S_r5)/2
row6=(r6+S_r6)/2
row7=(r7+S_r7)/2
row8=(r8+S_r8)/2
row9=(r9+S_r9)/2

[~,o1]=sort(row1,'descend');
[~,o2]=sort(row2,'descend');
[~,o3]=sort(row3,'descend');
[~,o4]=sort(row4,'ascend');
[~,o5]=sort(row5,'ascend');
[~,o6]=sort(row6,'ascend');
[~,o7]=sort(row7,'ascend');
[~,o8]=sort(row8,'ascend');
[~,o9]=sort(row9,'ascend');
