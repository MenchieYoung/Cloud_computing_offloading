%function he=hf(he_X)
function he=hf(maxx,maxy)
global W_sum D_sum;    % ant;
% he=cell(2,ant);
% for i=1:ant
% he{1,i}=sum(lastW_distri{1,i})-W_sum;
% he{2,i}=sum(lastD_distri{1,i})-D_sum;
 he(1,1)=sum(maxx)-W_sum;
 he(1,2)=sum(maxy)-D_sum;

end