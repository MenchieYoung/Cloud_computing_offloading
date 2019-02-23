% function [F]=C(x1,x2) %目标函数
% F=-(x1.^2+2*x2.^2-0.3*cos(3*pi*x1)-0.4*cos(4*pi*x2)+0.7);

%分配子载波和任务的数据量和计算量。 
%变量是子载波向量W，数据量向量D，假设计算量向量C=gama0*D= 设gama 令 C=(gama/alpha)*D
%目标函数Q=sum(Q(m))
%function [Q]= Q(W,D,alpha,m,a1,b1,a2,b2,a3,b3,N_0，gama0)    %alpha根据需要自己设置，alpha是时间代价的权重，beta是能耗代价 
                                %子任务数量m
                                %子任务发送功率范围[a1,b1]    单位mW
                                %节点对应信道增益的范围[a2,b2]
                                %节点计算能力的范围[a3,b3]
                                    %%怎么设置，有点纠结，毕竟会影响解的范围，不超过所有子任务计算量之和
                                %N_0:噪声功率   大约把SNR保持在20-50？
                                %C=gama0*D  
                                %SNR∈[a4,b4] 
                                
  function [Q]= Q(W,D,Mk,lamdak)   %Mk,lamdak应该是个向量，是随外迭代里次数而变化的 %Mk是罚因子


  % alpha=rand ;    %可以写到函数里
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%k可以移到主函数里吧？%%%%%%%%%%%%%%%%%%%

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%k可以移到主函数里吧？%%%%%%%%%%%%%%%%%%%
%    W0=diag(W);
%    D0=diag(D);
%    Q=(alpha+beta*P).*D0   log2(1+SNR).*W0  %w我太傻了，matlab有点乘点除啊
%    Q(m)=D(m)*(alpha+beta*P(m))/W(m)*log2(1+P(m)*g(m)/N_0)+alpha*gama0*D(m)/f(m)  
   %节点k计算能力f(k),假设一个子任务给一个节点(不然划分也没啥意义了)，那么f(k)就是f(m),反正F随机生成吧
%     for i=1:1:m
%     Q(i)=D(i)*(alpha+beta*P(i))/W(i)*log2(1+SNR(i))+alpha*gama0*D(i)/f(i)  %%感觉还是尽量用向量表示，不然太太太慢了
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global alpha P f SNR beta gama0 rho W_sum D_sum ;
    Q_obj=(D.*(alpha+beta*P))./(rho*(W.*log2(1+SNR)))+(alpha*gama0*D)./f;    
    %+ 1000*(D_sum-sum(D))+1000*(W_sum-sum(W)));
    Q=-(sum(Q_obj)+ Mk/2*((W_sum-sum(W))^2+(D_sum-sum(D))^2)- (lamdak*[W_sum-sum(W);D_sum-sum(D)]));   %因为要求sum（Q）+abs()+abs()的最小值，但蚁群算法求得是最大值,取个负号
                                       %感觉这里不对，我的h=W_sum-sum(maxx)
                                       %感觉还是要改hf，因为我的辅助函数是对所有蚂蚁都有要求的，固然h也是对所有蚂蚁而言
                                       
                               %abs（）最小值=0，也就是等式约束  
                               %真正的传输代价是sum（Q）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%     global f SNR gama0 rho W_sum D_sum ;
%     Q_obj=D./(rho*(W.*log2(1+SNR)))+(gama0*D)./f;  
%     Q=-(sum(Q_obj)+ Mk/2*((W_sum-sum(W))^2+(D_sum-sum(D))^2)- (lamdak*[W_sum-sum(W);D_sum-sum(D)]));