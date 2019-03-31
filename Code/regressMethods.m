  function [outCT] =regressMethods(iqa_L,gt_L)
  [s1, s2]=size(iqa_L);
  outCT=zeros(s1,s2);
  start=[0.0,0.1,0.0,0.0,0.0]';
  modelFun=@(b,x) b(1).*((1/2)-1./(1+exp(b(2).*(x-b(3)))))+b(4).*x+b(5);

  
for ii=1:s1 
nlmColDist=fitnlm(iqa_L(ii,:)',gt_L,modelFun,start);
outCT(ii,:)=predict(nlmColDist,iqa_L(ii,:)')';
end
  
  end
  
  