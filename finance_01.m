% Simulation 01
% author:Wei-Jian Huang uid:221081203010

rng('default')
mu = 0.05;          % 涨幅
sigma = 0.25;        % 波动率
S0 = 12;            % 起始价格
K = 15;             % 行使日价格
N = 10000;          % N条随机过程
M = 180;            % 观测点数
T = 0.5;              % 时间段
dt = T/M;           % 间隔时间
dW = randn(M,N);    % 标准正态分布
delta = [0:dt:0.5];  % 时间点
Seuler = zeros(M,N); % 欧拉法构造随机过程
Smilstein=zeros(M,N); % milstein构造随机过程
Seuler(1,:) = S0;       
Smilstein(1,:)=S0;
increment = 200; 

%% 构建随机过程
for i = 1:M
    %Euler Maruyama
    Seuler(i+1,:) = Seuler(i,:) + mu.*Seuler(i,:).*dt + sigma.*Seuler(i,:).*dW(i,:).*sqrt(dt); 
    %Milstein
    Smilstein(i+1,:) = Smilstein(i,:) + mu.*Smilstein(i,:).*dt - 0.5.*Smilstein(i,:)*((sigma)^2).*dt+ sigma.*Smilstein(i,:).*dW(i,:).*sqrt(dt) +0.5.*Smilstein(i,:).*((sigma)^2.*(dW(i,:).*sqrt(dt)).^2); %Milstein Scheme
end

%% 各取5个路径
figure;
subplot(2,1,1);
plot(delta, Seuler(:,1), 'r-',delta, Seuler(:,2), 'r-',delta, Seuler(:,3), 'r-',delta,Seuler(:,4),'r-',delta,Seuler(:,5), 'r-', delta,Seuler(:,6),'r-')
legend('Sample Paths for Euler M')
xlabel('Time')
ylabel('Stock Price')
title('欧拉法构造随机过程')

subplot(2,1,2);
plot(delta, Smilstein(:,1), 'y-',delta, Smilstein(:,2), 'y-',delta, Smilstein(:,3), 'y-',delta,Smilstein(:,4),'y-',delta,Smilstein(:,5), 'y-', delta,Smilstein(:,6),'y-')
legend('Sample Paths for Milstein')
xlabel('Time')
ylabel('Stock Price')
title('Milstein法构造随机过程')

%% 计算解析解和数值解的误差
Optioneuler = zeros(N/increment,1); %  期权价格的样本集合（euler）
Optionmilstein = zeros(N/increment,1); %  期权价格的样本集合（milstein）
Paths1 = zeros(N/increment,1); %  samples  of euler
Paths2 = zeros(N/increment,1); %   samples of milstein
for j = increment:increment:N
    eulerstock = Seuler(end,1:j)'; %  euler估计的股票价格
    milsteinstock = Smilstein(end,1:j)'; % milstein估计的股票价格
    Optioneuler(j/increment) = mean(max(eulerstock - K,0)); % 每j条估计路径的平均值
    Optionmilstein(j/increment) = mean(max(milsteinstock - K,0)); % 每j条估计路径的平均值
    Paths1(j/increment) = j; 
    Paths2(j/increment) = j;
end

% Black-Scholes的解析解
d1 = (log(S0/K) + (mu + 0.5*sigma^2)*T)/(sigma*sqrt(T)); 
d2 = d1 - sigma*sqrt(T);
N1 = 0.5*(1+erf(d1/sqrt(2)));
N2 = 0.5*(1+erf(d2/sqrt(2)));
C = S0*N1-K*exp(-mu*T)*N2;
g=length(Paths1);
C1=zeros(g+1,1);
C1(:,:)=C;

figure;
plot(Paths1,Optioneuler,Paths2,Optionmilstein,[0; Paths1],C1)
hold on
plot(0,C, '-o')
ylim([C-0.1 C+0.2]);
xlabel('Number of Sample Paths N')
ylabel('Call Option Price')
legend('Euler Maruyama Price','Milstein Price','Black Scholes Price')
title('Comparison of Option Price between Euler-M,Milstein and BS')

