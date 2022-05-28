clear all;
close all;
clc;

M = 10;  
N = 10;
c = 3e8;
f = 22e9;
lamda = c/f;             %波长
d = lamda/2;             %阵元间距
theta0 = 0;             %入射方向
fai0 = 90;
u0 = sind(theta0)*cosd(fai0);
v0 = sind(theta0)*sind(fai0);
bujing = 0.001;
range = -1 :bujing :1;       %角度
bs = 1- range(1)/bujing;
Mk = round(1*M*N);                 %开启的阵元数
str1 = ones(1,Mk);
str2 = zeros(1,M*N - Mk);
Sta = [str1 str2];
Sta = Sta(randperm(M*N));
Sta(1) = 1;               %阵列孔径约束
Sta(M) = 1;
Sta(M*(N-1)+1) = 1;
Sta(M*N) = 1;
  
N_ones = sum(Sta);
if(N_ones > Mk)           %阵列稀疏率约束
    nn1 = find(Sta==1);
    nn11 = nn1((nn1~=1)&(nn1~=M)&(nn1~= M*(N-1)+1)&(nn1~=M*N));
    nn111 = nn11(randperm(numel(nn11),(N_ones - Mk)));
    for rr = 1 :(N_ones - Mk )
        Sta(nn111(rr)) = 0;
    end
elseif(N_ones < Mk)
    nn22 = find(Sta==0);
    nn222 = nn22(randperm(numel(nn22),(Mk - N_ones)));
    for rr = 1 :(Mk - N_ones)
        Sta(nn222(rr)) = 1;
    end
end

% Sta = gggg;
Stb = reshape(Sta,M,N);
fitness = zeros(length(range),length(range));

m=(0:M-1)'; % m行
n=(0:N-1);  % n列

tao_W=(m*d*u0+n*d*v0)/c; 
W = Stb.*exp(1j*2*pi*f*tao_W);
W1=reshape(W,[],1);     

for v = range
    for u = range
        tao = (m*d*u+n*d*v)/c; 
        V = exp(1j*2*pi*f*tao);     % 方向矢量  
        V1=reshape(V,[],1);      
        fitness(round((u/bujing) + bs),round((v/bujing) + bs)) = W1'*V1;
    end
end

fitness1 = abs(fitness);
fitness = 20*log10(abs(fitness./max(max(fitness))));
figure;
mesh(range,range,fitness);
set(gca,'FontSize',15);
axis([-inf,inf,-inf,inf,-30 0])
xlabel('v','fontsize',15);
ylabel('u','fontsize',15);
zlabel('(幅度方向图（dB)','fontsize',15);
figure;
mesh(range,range,fitness1);
xlabel('v');
ylabel('u');
zlabel('(幅度方向图)');
figure;
plot(range,fitness(:,round((v0/bujing) + bs)));
xlabel('u');
ylabel('(幅度（dB)');

figure;
plot(range,fitness(round((u0/bujing) + bs),:));
xlabel('v');
ylabel('(幅度（dB)');

max1 = -50;
max2 = -50;
max3 = -50;
max4 = -50;
result1 = -50;
result2 = -50;

if(u0 < 0)
    for bb1 = u0 : bujing : 1-bujing
      if (fitness(round((bb1/bujing) + bs),round((v0/bujing) + bs))<= fitness(round((bb1/bujing) + bs)+1,round((v0/bujing) + bs)) )
            range2 = bb1; 
            if(2*u0-range2 >= -1)
                max1 = max(max(fitness(1 :round(((2*u0-range2)/bujing) + bs) ,:)));
                max2 = max(max(fitness(round((range2/bujing) + bs):end ,:)));
                result1  = max(max1,max2);
            else
                result1  = max(max(fitness(round((range2/bujing) + bs) :round(((2*u0-range2+2)/bujing)+ bs) ,:)));
            end
            break;   
      end
    end
else
    for bb2 = u0 : -bujing : -1+bujing
      if (fitness(round((bb2/bujing) + bs),round((v0/bujing) + bs))<= fitness(round((bb2/bujing) + bs)-1,round((v0/bujing) + bs)) )
            range3 = bb2; 
            if(2*u0-range3 <= 1)
                max1 = max(max(fitness(round(((2*u0-range3)/bujing) + bs):end ,:)));
                max2 = max(max(fitness(1:round((range3/bujing) + bs) ,:)));
                result1  = max(max1,max2);
            else
                result1  = max(max(fitness(round(((2*u0-range3-2)/bujing)+ bs):round((range3/bujing) + bs) ,:)));
            end
            break;   
      end
    end       
end


if(v0 < 0)
    for bb3 = v0 : bujing : 1-bujing
      if (fitness(round((u0/bujing) + bs),round((bb3/bujing) + bs))<= fitness(round((u0/bujing) + bs),round((bb3/bujing) + bs)+1) )
            range4 = bb3; 
            if(2*v0-range4 >= -1)
                max3 = max(max(fitness(:,1 :round(((2*v0-range4)/bujing) + bs))));
                max4 = max(max(fitness(:,round((range4/bujing) + bs):end )));
                result1  = max(max4,max3);
            else
                result1  = max(max(fitness(:,round((range4/bujing) + bs) :round(((2*v0-range4+2)/bujing)+ bs) )));
            end
            break;   
      end
    end
else
    for bb4 = v0 : -bujing : -1+bujing
      if (fitness(round((u0/bujing) + bs),round((bb4/bujing) + bs))<= fitness(round((u0/bujing) + bs),round((bb4/bujing) + bs)-1) )
            range5 = bb4; 
            if(2*v0-range5 <= 1)
                max3 = max(max(fitness(:,round(((2*v0-range5)/bujing) + bs):end)));
                max4 = max(max(fitness(:,1:round((range5/bujing) + bs))));
                result2  = max(max3,max4);
            else
                result2  = max(max(fitness(:,round(((2*v0-range5-2)/bujing)+ bs):round((range5/bujing) + bs))));
            end
            break;   
      end
    end       
end
result = max(result1 , result2);
%  