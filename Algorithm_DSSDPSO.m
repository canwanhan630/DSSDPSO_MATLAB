clear all;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ʼ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NP = 150;                     %Ⱥ����������
M = 20;                       %��������
N = 20;                       %��������
SPrate = 0.5;                 %ϡ����
Mk = round(SPrate*M*N);       %������Ԫ��
T = 500;                      %��������
c1 = 1.5;                     %ѧϰ����1
w = 0.8;                      %����Ȩ��
vmax = 10;                    %�ٶ����ֵ
vmin = -10;                   %�ٶ���Сֵ


dimQ = round(M*N/2);           %ȫ�ֲ���ѡ��ά����
sol = 3;                       %ȫ�ֲ���ѡ���������


thr = 0.9;                     %�ֲ������������         
WINpara = 4/3;                 %���Ų���
LOSEpara = 3/4;                %��������

Tpb = 5;                      %�ֲ����ų�������
Txb = 20;                     %�����ռ�������
MUT_pro = 1/5;                % �������



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ʼ����Ⱥ����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = zeros(NP,M*N);          %������ʼ��Ⱥ
for lo = 1 : NP
  str1 = ones(1,Mk);
  str2 = zeros(1,M*N - Mk);
  Sta = [str1 str2];
  Sta = Sta(randperm(M*N));
  
  Sta(1) = 1;               %���п׾�Լ��
  Sta(M) = 1;
  Sta(M*(N-1)+1) = 1;
  Sta(M*N) = 1;
  
  N_ones = sum(Sta);
  if(N_ones > Mk)           %����ϡ����Լ��
      temp11 = find(Sta==1);
      temp12 = temp11((temp11~=1)&(temp11~=M)&(temp11~= M*(N-1)+1)&(temp11~=M*N));
      temp13 = temp12(randperm(numel(temp12),(N_ones - Mk)));
      for rr = 1 :(N_ones - Mk )
          Sta(temp13(rr)) = 0;
      end
  elseif(N_ones < Mk)
      temp21 = find(Sta==0);
      temp22 = temp21(randperm(numel(temp21),(Mk - N_ones)));
      for rr = 1 :(Mk - N_ones)
          Sta(temp22(rr)) = 1;
      end
  end        
  x(lo,:) = Sta;     
end

v = zeros(NP,M*N);
for lo = 1 : NP
  Stv = rand(1,M*N)*(vmax - vmin) + vmin; 
  v(lo,:) = Stv;     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ʼ����������λ��������ֵ%%%%%%%%%%%%%%%%%%%%%%%%

p = x;
pbest = ones(NP,1);
for i = 1 :NP
    AAAAWD = func1(x(i,:));
    pbest(i) = AAAAWD;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ʼ��ȫ������λ��������ֵ%%%%%%%%%%%%%%%%%%%%%%%%

g = ones(1,M*N);
gbest = inf;
for i = 1 :NP
    if(pbest(i) < gbest)
        g = p(i,:);
        gbest = pbest(i);
    end
end
gb = ones(1,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : T
    
    sTag = zeros(1,NP);
    rTag = zeros(1,NP);
    for j = 1 : NP       

%���¸�������λ�ú�����ֵ-----------------------------
        xfun1 = func1(x(j,:));
        if( xfun1 < pbest(j) )
            p(j,:) = x(j,:);
            pbest(j) = xfun1;
            sTag(j) = 0;
        else
            sTag(j) = sTag(j) + 1;            
        end
        
%����ȫ������λ�ú�����ֵ-----------------------------      
 
        if( pbest(j) < gbest )
            g = p(j,:);
            gbest = pbest(j);
        end       
        
%����λ�ú��ٶ�ֵ------------------------------------  

        ZZQ = linspace(1,NP,NP);                                      %�ۺ�����ѧϰ����
        solu_PG = ZZQ(randperm(numel(ZZQ),sol));
        if(j == 1)
            solu_PG_Plus = [solu_PG NP 2];
        elseif(j == NP)
            solu_PG_Plus = [solu_PG NP-1 1];
        else
            solu_PG_Plus = [solu_PG j-1 j+1];
        end
        solu_PG_Plus_Value = zeros(1,length(solu_PG_Plus));
        for yyx = 1 :length(solu_PG_Plus)
            solu_PG_Plus_Value(yyx) = pbest(solu_PG_Plus(yyx));
        end
        solu_PG_Find1 = find(solu_PG_Plus_Value==min(solu_PG_Plus_Value));
        solu_PG_Find2 = solu_PG_Plus(solu_PG_Find1(1));    

        ZWD = linspace(1,M*N,M*N);
        XXWD = ZWD(randperm(numel(ZWD),dimQ));
        xxrand = rand;
        v(j,:) = w*v(j,:) + c1*xxrand*( p(j,:) - x(j,:));
        for wd = 1 : dimQ
            v(j,XXWD(wd)) = w*v(j,XXWD(wd)) + c1*xxrand*( p(solu_PG_Find2,XXWD(wd)) - x(j,XXWD(wd)));
        end
        ss = zeros(1,M*N);
        ff = x(j,:);
        for k = 1 : M*N
            ss(k) = 1/(1 + exp( -v(j,k) ));
            r = rand;
            if r < ss(k)
                x(j,k) = 1;
            else
                x(j,k) = 0;
            end
        end
        
        
        
        if(i >= thr*T)                                 %�ֲ������������
            OPT_pro = 1/12;                            %��ʼ�ֲ������������
            while OPT_pro >= 0.05
                temp_x = zeros(1,M*N);
                for k = 1 : M*N                        
                    r = rand;
                    if r < OPT_pro
                        temp_x(k) = ~x(j,k);
                    else
                        temp_x(k) = x(j,k);
                    end
                end
                if( func1(temp_x) < func1(x(j,:)) )
                    x(j,:) = temp_x;
                    OPT_pro = WINpara*OPT_pro;
                else
                    OPT_pro = LOSEpara*OPT_pro;
                end
            end
        end
        
        
        if((sTag(j) > Tpb )&&(rTag(j) > Txb ))          %����������
            for k = 1 : M*N
                r = rand;
                if r < MUT_pro
                    x(j,k) = ~x(j,k);
                else
                    x(j,k) = x(j,k);
                end
            end   
        end
                     
        x(j,1) = 1;               %���п׾�Լ��
        x(j,M) = 1;
        x(j,M*(N-1)+1) = 1;
        x(j,M*N) = 1;
        
        NN_ones = sum(x(j,:));
        if(NN_ones > Mk)           %����ϡ����Լ��
            temp31 = find(x(j,:)==1);
            temp32 = temp31((temp31~=1)&(temp31~=M)&(temp31~= M*(N-1)+1)&(temp31~=M*N));
            temp33 = temp32(randperm(numel(temp32),(NN_ones - Mk)));
            for rr = 1 :(NN_ones - Mk )
                x(j,temp33(rr)) = 0;
            end
        elseif(NN_ones < Mk)
            temp41 = find(x(j,:)==0);
            temp42 = temp41(randperm(numel(temp41),(Mk - NN_ones)));
            for rr = 1 :(Mk - NN_ones)
                x(j,temp42(rr)) = 1;
            end
        end
        ff2 = ff - x(j,:);
        if(sum(ff2) < round(M*N*MUT_pro))
            rTag(j) = rTag(j) + 1 ;
        else
            rTag(j) = 0;
        end
        
        
 %�߽���������--------------------------------------
 
        for ii = 1 : M*N
            if ( v(j,ii) > vmax  ||   v(j,ii) < vmin   )
                v(j,ii) = rand*(vmax - vmin) + vmin;
            end
        end
    end
        
%��¼����ȫ������ֵ----------------------------------
        
    gb(i) = gbest;
end
figure;
plot(gb);
xlabel('��������')
ylabel('��Ӧ��ֵ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��Ӧ�Ⱥ�����PSLL��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        
function result = func1(X)

M = 20;  
N = 20;
c = 3e8;
f = 22e9;
lamda = c/f;             %���䲨��
d = lamda/2;             %��Ԫ���
theta0 = 0;              %���䷽��
fai0 = 90;
u0 = sind(theta0)*cosd(fai0);
v0 = sind(theta0)*sind(fai0);
bujing = 0.01;
range = -1 :bujing :1;       %�Ƕ�
bs = 1- range(1)/bujing;
Stb = reshape(X,M,N);
fitness = zeros(length(range),length(range));

m=(0:M-1)'; % m��
n=(0:N-1);  % n��

tao_W=(m*d*u0+n*d*v0)/c; 
W = Stb.*exp(1j*2*pi*f*tao_W);
W1=reshape(W,[],1);     

for v = range
    for u = range
        tao = (m*d*u+n*d*v)/c; 
        V = exp(1j*2*pi*f*tao);     % ����ʸ��  
        V1=reshape(V,[],1);      
        fitness(round((u/bujing) + bs),round((v/bujing) + bs)) = W1'*V1;
    end
end

fitness = 20*log10(abs(fitness./max(max(fitness))));
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

end        