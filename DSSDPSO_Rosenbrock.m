clear all;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%初始化%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NP = 200;     %群体粒子数量
M = 1;
N = 10;
T = 1000;      %迭代次数
c1 = 1.5;     %学习因子1
w = 0.8;      %惯性权重
vmax = 10;
vmin = -10;
xmax = 30;
xmin = -30;
wei = round(M*N/2);   
yx = 3;

afa =  0.9;
wingama = 4/3;
losegama = 3/4;

tt = 5;       %局部最优持续代数
nn = 20;      %搜索空间计算代数
gdfg = 1/5;  % 变异概率

cishu = 50;
AVERANGE  = zeros(1,cishu);

for kk = 1 : cishu
    jilu2 = kk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%初始化种群个体%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = zeros(NP,M*N);      %产生初始种群
for lo = 1 : NP
   Sta = xmin + (xmax - xmin)*rand(1,M*N);  
  x(lo,:) = Sta;     
end

v = zeros(NP,M*N);
for lo = 1 : NP
  Stv = rand(1,M*N)*(vmax - vmin) + vmin; 
  v(lo,:) = Stv;     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%初始化个体最优位置与最优值%%%%%%%%%%%%%%%%%%%%%%%%

p = x;
pbest = ones(NP,1);
for i = 1 :NP
    AAAAWD = func1(x(i,:));
    pbest(i) = AAAAWD;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%初始化全局最优位置与最优值%%%%%%%%%%%%%%%%%%%%%%%%

g = ones(1,M*N);
gbest = inf;
for i = 1 :NP
    if(pbest(i) < gbest)
        g = p(i,:);
        gbest = pbest(i);
    end
end
gb = ones(1,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%迭代%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : T
    jilu = i;
    sTag = zeros(1,NP);
    rTag = zeros(1,NP);
    for j = 1 : NP       

%更新个体最优位置和最优值-----------------------------
        xfun1 = func1(x(j,:));
        if( xfun1 < pbest(j) )
            p(j,:) = x(j,:);
            pbest(j) = xfun1;
            sTag(j) = 0;
        else
            sTag(j) = sTag(j) + 1;            
        end
        
%更新全局最优位置和最优值-----------------------------      
 
        if( pbest(j) < gbest )
            g = p(j,:);
            gbest = pbest(j);
        end       
        
%更新位置和速度值------------------------------------         
        ZZQ = linspace(1,NP,NP);                                      %综合学习策略
        YouXiu = ZZQ(randperm(numel(ZZQ),yx));
        if(j == 1)
            YouXiuPlus = [YouXiu NP 2];
        elseif(j == NP)
            YouXiuPlus = [YouXiu NP-1 1];
        else
            YouXiuPlus = [YouXiu j-1 j+1];
        end
        YouXiuPlusValue = zeros(1,length(YouXiuPlus));
        for yyx = 1 :length(YouXiuPlus)
            YouXiuPlusValue(yyx) = pbest(YouXiuPlus(yyx));
        end
        YouXiuFind = find(YouXiuPlusValue==min(YouXiuPlusValue));
        YouXiuFind2 = YouXiuPlus(YouXiuFind(1));    

        ZWD = linspace(1,M*N,M*N);
        XXWD = ZWD(randperm(numel(ZWD),wei));
        xxrand = rand;
        v(j,:) = w*v(j,:) + c1*xxrand*( p(j,:) - x(j,:));
        for wd = 1 : wei
            v(j,XXWD(wd)) = w*v(j,XXWD(wd)) + c1*xxrand*( p(YouXiuFind2,XXWD(wd)) - x(j,XXWD(wd)));
        end
        ff = x(j,:);
        x(j,:) = x(j,:) + v(j,:);
    
        if(i >= afa*T)
            JUBUgdfg = 1/12;  %初始局部搜索变异概率
            while JUBUgdfg >= 0.05
                jubux = zeros(1,M*N);
                for k = 1 : M*N                        %局部搜索变异策略
                    r = rand;
                    if r < JUBUgdfg
                        jubux(k) = x(j,k) + 0.05*(xmin + (xmax - xmin)*rand);
                        if ( jubux(k) > xmax  ||   jubux(k) < xmin   )
                            jubux(k) = x(j,k);
                        end
                    else
                        jubux(k) = x(j,k);
                    end
                end
                if( func1(jubux) < func1(x(j,:)) )
                    x(j,:) = jubux;
                    JUBUgdfg = wingama*JUBUgdfg;
                else
                    JUBUgdfg = losegama*JUBUgdfg;
                end
            end
        end
        
        
        if((sTag(j) > tt )&&(rTag(j) > nn ))   %个体变异策略
            for k = 1 : M*N
                r = rand;
                if r < gdfg
                    x(j,k) = x(j,k) + 0.3*(xmin + (xmax - xmin)*rand);
                    if ( x(j,k) > xmax  ||   x(j,k) < xmin   )
                        x(j,k) = x(j,k);
                    end
                else
                    x(j,k) = x(j,k);
                end
            end   
        end
                     
        ff2 = ff - x(j,:);
        if(sum(ff2) < round(M*N*gdfg))
            rTag(j) = rTag(j) + 1 ;
        else
            rTag(j) = 0;
        end
        
        
 %边界条件处理--------------------------------------
 
        for ii = 1 : M*N
            if ( v(j,ii) > vmax  ||   v(j,ii) < vmin   )
                v(j,ii) = rand*(vmax - vmin) + vmin;
            end
        end
        for ii = 1 : M*N
            if ( x(j,ii) > xmax  ||   x(j,ii) < xmin   )
                x(j,ii) = rand*(xmax - xmin) + xmin;
            end
        end
    end
        
%记录历代全局最优值----------------------------------
        
    gb(i) = gbest;
end
 
AVERANGE(kk) = gb(T);
end

AVERANGE_mean = mean(AVERANGE);
AVERANGE_var = var(AVERANGE);
AVERANGE_min = min(AVERANGE);
% figure;
% plot(gb);
% xlabel('迭代次数')
% ylabel('适应度值')

function [y] = func1(xx)

d = length(xx);
sum = 0;
for ii = 1:(d-1)
	xi = xx(ii);
	xnext = xx(ii+1);
	new = 100*(xnext-xi^2)^2 + (xi-1)^2;
	sum = sum + new;
end

y = sum;

end
