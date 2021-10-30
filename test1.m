clear all;
clc
data_bits_in=rand(2,6)<0.5;
[m,n]=size(data_bits_in);
tx_64QAM=Mapping_64QAM(data_bits_in,m,n)
data_bits_out=Demapping_64QAM(tx_64QAM,m,n)
data_bits_in-data_bits_out




function [tx_64QAM]=  Mapping_64QAM (tx_bits,m,n)
num=m*n;
map1=[-7,-7,-7,-7,-7,-7,-7,-7,    -5,-5,-5,-5,-5,-5,-5,-5,   -1,-1,-1,-1,-1,-1,-1,-1, -3,-3,-3,-3,-3,-3,-3,-3,  7,7,7, 7, 7,7,7,7,     5 ,5,5 ,5 , 5 ,5,5 ,5,     1,1,1,1,1,1,1,1,    3,3,3,3,3,3,3,3];
map2=[-7,-5,-1,-3,7,5,1,3,     -7, -5, -1 ,-3, 7, 5, 1, 3,     -7, -5, -1, -3, 7, 5, 1, 3,     -7, -5, -1, -3, 7, 5, 1, 3,      -7, -5, -1, -3, 7, 5 ,1, 3,      -7, -5, -1, -3, 7 ,5, 1, 3 ,     -7 ,-5, -1, -3, 7 ,5 ,1 ,3 ,    -7, -5, -1, -3, 7, 5 ,1 ,3];

%初始化
Dec_ac=zeros(num/6,1);
ac=zeros(num/6,1);
as=zeros(num/6,1);
tx_64QAM=zeros(num/6,1);

 
 for N=1:num/6
    Dec_ac(N,:)=32*tx_bits(N,1)+16*tx_bits(N,2)+8*tx_bits(N,3)+4*tx_bits(N,4)+2*tx_bits(N,5)+tx_bits(N,6);    %第N行求十进制 
    %同相分量ac
    ac(N)=map1(Dec_ac(N)+1);%MATLAB索引从1开始
    %正交分量as
    as(N)=map2(Dec_ac(N)+1);%MATLAB索引从1开始  
    tx_64QAM(N,1)=ac(N)+1i*as(N);    %返回一个虚数
 end
end



function [bin] = Demapping_64QAM(rx_64QAM,m,n)
num=m*n;
yc=real(rx_64QAM);% 取实部
ys=imag(rx_64QAM);%取虚部
bin=ones(num/6,6);

map1=[-7,-5,-1,-3,7,5,1,3];
% 通过map1把坐标映射到demapping中的二进制数
demapping=[
    0 0 0;
    0 0 1;
    0 1 0;
    0 1 1;
    1 0 0;
    1 0 1;
    1 1 0;
    1 1 1];
for N=1:num/6
%最大似然门限判决 
    if(yc(N)<=-6)
        yc(N)=-7;     
    elseif(yc(N)>-6&&yc(N)<=-4)
        yc(N)=-5;
    elseif(yc(N)>-4&&yc(N)<=-2)
        yc(N)=-3;
    elseif(yc(N)>-2&&yc(N)<=0)
        yc(N)=-1; 
    elseif(yc(N)>0&&yc(N)<=2)
        yc(N)=1;
    elseif(yc(N)>2&&yc(N)<=4)
        yc(N)=3;
    elseif(yc(N)>4&&yc(N)<=6)
        yc(N)=5;
    elseif(yc(N)>6)
        yc(N)=7;
    end
 
 
    if(ys(N)<=-6)
        ys(N)=-7;
    elseif(ys(N)>-6&&ys(N)<=-4)
        ys(N)=-5;
    elseif(ys(N)>-4&&ys(N)<=-2)
        ys(N)=-3;
    elseif(ys(N)>-2&&ys(N)<=0)
        ys(N)=-1;
    elseif(ys(N)>0&&ys(N)<=2)
        ys(N)=1;
    elseif(ys(N)>2&&ys(N)<=4)
        ys(N)=3;
    elseif(ys(N)>4&&ys(N)<=6)
        ys(N)=5;
    elseif(ys(N)>6)
        ys(N)=7;
    end
end  %对应for循环
%     解映射
    index=cell(1,N);
    index(yc==-7)={[0 0 0]};
    index(yc==-5)={[0 0 1]};
    index(yc==-1)={[0 1 0]};
    index(yc==-3)={[0 1 1]};
    index(yc==7)={[1 0 0]};
    index(yc==5)={[1 0 1]};
    index(yc==1)={[1 1 0]};
    index(yc==3)={[1 1 1]};
    a1=reshape(cell2mat(index),[3,N])';
    
    index=cell(1,6);
    index(ys==-7)={[0 0 0]};
    index(ys==-5)={[0 0 1]};
    index(ys==-1)={[0 1 0]};
    index(ys==-3)={[0 1 1]};
    index(ys==7)={[1 0 0]};
    index(ys==5)={[1 0 1]};
    index(ys==1)={[1 1 0]};
    index(ys==3)={[1 1 1]};
    a2=reshape(cell2mat(index),[3,N])';
    bin=[a1,a2];
end