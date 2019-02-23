% 蚁群算法求函数最大值的程序
%function max_need=antcolony_offload(m,)  
function [lastW_distri,lastD_distri,maxx,maxy,minvalue]=antcolony_offload(m,X,Mk,lamdak,k,times)
%global W_sum;
%global alpha P f SNR beta  gama0  W_sum D_sum rho;
global ant;
%ant=200;   % 蚂蚁数量
%times=50; % 蚂蚁移动次数
rou=0.8; % 信息素挥发系数
p0=0.2; % 转移概率常数
lower_1=zeros(1,m); % 设置搜索范围  &横   %W向量
upper_1=20.*ones(1,m); %   横
lower_2=zeros(1,m); %  纵     %D
upper_2=30.*ones(1,m); %   纵
%subcarrier=zeros(1,m);
% X{1,:}=zeros(1,m);   %有错 还没管  %其实X{1,:}这样是X{1,1},是个向量
% X{2,:}=zeros(1,m);
% X=cell(2,ant);
% %如果定义成矩阵，那么想让元素为cell且提取出这个cell的内容用X{}，如果定义成cell，那么提取cell用（），提取cell内容为矩阵用{}
p1=zeros(times,ant);
tau_best=zeros(1,times);
minvaluet=zeros(1,times);
tau=zeros(1,ant);
%N=512;    %子载波个数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%调试用%%%%%%%%%%%
  %a1 b1 a3 b3 a4 b4 gama0 W_sum D_sum Q beta P f SNR; 
%     m=4;    %从第六个点开始突变。。。
%     alpha=0.5;
%     a1=50;     %mW
%     b1=100;    %mW
%     a3=8;  %8000kB/s  % 8000khz  %8Mhz
%     b3=10;  %Mhz
%     a4=20;
%     b4=35;
%     gama0=20;
%     W_sum=20;       %总子载波数 %先考虑等式，学姐论文开始也一直认为是=   %拉格朗日（回去找之前论文，没找到）  罚函数(万方论文)           %总带宽吧
%     D_sum=30;        
% % Q0=zeros(1,m);  %初始化
%     beta=1-alpha ;
%     rho=2^(-1/2);
%     P = a1 + (b1-a1).*rand([1 m]);  %在区间（a，b）内产生均匀分布的随机数   %m行 n列
% % %    g = a2 + (b2-a2).*rand([1 m]);
%     f = a3 + (b3-a3).*rand([1 m]);
%     SNR = a4 + (b4-a4).*rand([1 m]);
%    Mk=100;
%    lamdak=5;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:ant 
%   X{1,i}=(lower_1+(upper_1(1,1)-lower_1(1,1)).*rand(1,m)); %随机设置每只蚂蚁的初值位置  %横 W
%   X{2,i}=(lower_2+(upper_2(1,1)-lower_2(1,1)).*rand(1,m)); %纵  D
%   
%   tau(i)= Q(X{1,i},X{2,i},Mk,lamdak); %第i只蚂蚁的信息量   
%                            %%所在位置函数值越大，信息量越大
% end
for i=1:ant
    tau(i)=Q(X{1,i},X{2,i},Mk,lamdak);  %lamdak(:,i)列向量
end

for t=1:times  % 第t次移动
    lamda=1/t; %步长系数，随移动次数增大而减少  %%就像二值逼近，分的越来越细，因为越来越接近所需值  %%%可以考虑改变此关系函数？
    [tau_best(t),bestindex]=max(tau); %第t次移动的最优值及其位置  %所有蚂蚁所在位置F最大处  %第一次移动，t=1时，也就是初始位置
    %taubestt(t)=tau(bestindex);  %和taubest一样，另外它也没错啊
    minvaluet(t)=-Q(X{1,bestindex},X{2,bestindex},Mk,lamdak);
    %%%%%%%%%%%为什么不是相反数？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？
    %%%%%%%%这个bestindex是对i而言的没错，就该画这个图而不是tau_best的，tau参量是i，所以问题不在这
    %%%%%%%%%傻死了，tau是信息素，还有之前的印记，-minvaluet完全就是当前的目标函数值，当然不一样。只有第一次是相反数
    
    for i=1:ant    %第i只蚂蚁
      p1(t,i)=(tau(bestindex)-tau(i))/tau(bestindex); %最优值与第i只蚂蚁的值的差比
                                         %计算状态转移概率   %%差的越大，转移概率越大？
    end 
    
    for i=1:ant 
       if p1(t,i)<p0   %局部搜索  小于转移概率常数  %离最大值不远，不怎么要转移，就算lamda=1，也只走了不到1的步长
                                 
         temp1=X{1,i}+(2.*rand(1,m)-1)*lamda; %移动距离  %移动次数越多，步长越小（逼近法）  %rand应该是处于0-1，所以2*rand-1考虑了二维上移动的随机性
                                        %n维应该也一样，每一维都有正负两个方向 %有n个temp，要弄个矩阵
         temp2=X{2,i}+(2.*rand(1,m)-1)*lamda; 
       else 
      %全局搜索
         temp1=X{1,i}+(upper_1(1,1)-lower_1(1,1)).*(rand(1,m)-0.5);   %%还是会超出范围吧？比如X在上半部分，rand=1，会跑到上界的上面   后面有越界处理 超过的算边界
                                                      %有点像重新随机分配位置，但其实是在原来基础上随意重新无规律走..
         temp2=X{2,i}+(upper_2(1,1)-lower_2(1,1)).*(rand(1,m)-0.5); 
       end  
%%%%%%%%%%%%%%%%%%%%%% 越界处理%%%%%%%%%%%%%%%%%%%%%
       for j=1:1:m
         if temp1(j)<lower_1(j);    %其实lower里元素都一样，目前都是0
            temp1(j)=lower_1(j); 
         end 
         if temp1(j)>upper_1(j) 
            temp1(j)=upper_1(j); 
         end 
         if temp2(j)<lower_2(j) 
            temp2(j)=lower_2(j); 
         end 
         if temp2(j)>upper_2(j) 
            temp2(j)=upper_2(j); 
         end 
       end
%%%%%%%%%%%%%%%%%%%%%%% 
       if Q(temp1,temp2,Mk,lamdak)>Q(X{1,i},X{2,i},Mk,lamdak) % 判断蚂蚁是否移动 %保留那些走成功的，在此基础上随机乱走。没成功的再上一次基础上再次随机乱走，
         X{1,i}=temp1;         %%终于更新了蚂蚁位置
         X{2,i}=temp2; 
       end 
    end 
    
    for i=1:ant
       tau(i)=(1-rou).*tau(i)+Q(X{1,i},X{2,i},Mk,lamdak);  %更新信息量   %因为一些蚂蚁点更新了，F变了 
                                               %不管更没更新，所有蚂蚁都是加强了信息素，但是正确方向的蚂蚁Q更大，加强更明显
       
    end  
    
end 

 lastW_distri=X(1,:);   %ant个W向量（每个是个cell） %试了是用cell 也就是说用（）。但是我要的应该是内容所以应该是{}吧，但是{1，：}取出的是第一个
 lastD_distri=X(2,:);   %所以mulphr里X0不等于{lastW_distri;lastD_distri},因为人家已经是cell了= =
 t=(1:1:times);
%  figure,plot(t,tau_best(t)),title('信息素每次迭代最值(所有蚂蚁中)变化历程');hold on;
%  text(150,tau_best(150),num2str(k));
 plot(t,minvaluet(t)),title('辅助函数内迭代最值(所有蚂蚁中)变化历程');hold on;  %应该和上面那个是相反数才对
%text(10+2*k,minvaluet(10+2*k),num2str(k));
%text(k+1,minvaluet(k+1),num2str(k));
[~,max_index]=max(tau);       %好奇怪，找到信息素最强的对应的位置。。。为什么不直接找此时eval(f)里最大的？  %是一样的
maxx=X{1,max_index};              %W
maxy=X{2,max_index};              %D
maxvalue=Q(X{1,max_index},X{2,max_index},Mk,lamdak);   %和max_value有区别吗 %没有吧
minvalue=-maxvalue;   
%max_need=Q_need(maxx,maxy); 

% accu=W_sum/N; 
% for i=1:m
%   if mod(maxx(i),accu)>=accu/2
%     maxx(i)=ceil(maxx(i)/accu)*accu;
%   else
%     maxx(i)=floor(maxx(i)/accu)*accu;
%   end
% end
% subcarrier=maxx./accu;

%%end 