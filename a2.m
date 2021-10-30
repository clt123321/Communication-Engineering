clear all;
clc
%
%17206018 ChengLitao
%
% test Hamming code performance
u=[1 1 0 0];
P=[1 1 0; 0 1 1; 1 1 1 ;1 0 1];
I3=eye(3);
I4=eye(4);
G=[P,I4];
c=mod(u*G,2);
%c1=[1    1    0    1     0     0     0]
H=[I3,P'];
s=mod(c*H',2);
%s1=mod(c1*H',2)

%1.plot show the simulated and theoretical bit error rate (SER) versusE b /N 0 in dB.
%2.plot show the simulated and theoretical bit error rate (BER) versusE b /N 0 in dB.
%3.find cross point which the Hamming codeoffers improved performance over an uncoded system(BER)
% 4.estimate the value of E s /N 0 above which the system BER lies below 10^4
groupNum=1;
symbolNum=4*groupNum;
d=1;
Es=2.5*d*d;%enger of signal power(d^2)

sx=[-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5];
sy=[1.5,1.5,1.5,1.5,0.5,0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5,-1.5,-1.5,-1.5,-1.5];
SER=[];
BER=[];
%sigma=0.5;
for sigma=0.11:0.01:0.5
    %1.CreatbitSignal
    data_bits_in=rand(symbolNum,4)<0.5
    %2.Hamming encode
    code_bits_in=mod(data_bits_in*G,2)
    bit_in=(reshape(code_bits_in,[4,7*groupNum]))'


    %3.16-QAM modulation
    BinToDec=bit_in(:,1)*8+bit_in(:,2)*4+bit_in(:,3)*2+bit_in(:,4);
    symbol_in=MapToSymbol(BinToDec,7*groupNum);
    v=modulation(symbol_in,sx,sy,7*groupNum);

    %4.add AWGN
    noise=sigma.*randn(7*groupNum,2);
    v=v+noise;

    %5.16-QAM demodulation £¨optimum£©
    s=[sx' sy'];
    ds=zeros(7*groupNum,16);
    for i=1:16
        ds(:,i)=(v(:,1)-sx(i)).^2+(v(:,2)-sy(i)).^2;
    end
    [min_d,symbol_out]=min(ds,[],2);
    symbol_out;



    %6.Hamming decode(3 steps)
    bits_out=MapToBits(symbol_out,7*groupNum)
    code_bits_out=(reshape(bits_out,[7,symbolNum]))'

    
    

    %7calculate SER & BER
    a= (symbol_in'==symbol_out);
    SER=[SER,1-(sum(a(:)==1)/(7*groupNum))]

    b=bit_in==bits_out;
    BER=[BER,1-(sum(b(:)==1)/(4*7*groupNum))]
    
    
end

sigma=0.11:0.01:0.5;
N0=2*(sigma.^2);% En=noise RMS value^2=N0/2
theory_value=3*qfunc(d./(sqrt(2*N0)))-2.25.*(qfunc(d./(sqrt(2.*N0)))).^2;

%plot SER to Es/N0
x=10.*log10(Es./N0);
figure(1);
semilogy(x,SER,'-b.','MarkerSize',10);
hold on;
semilogy(x,theory_value,'-r.','MarkerSize',5);
legend("uncoding","theory value of uncoding");
title('SymbolErrorRate to Es/N0 in AWGN');
xlabel('Es/N0(dB)');
ylabel('SymbolErrorRate(SER)');
hold off;


%plot BER to Eb/N0
figure(2);
semilogy(10.*log10(Es./(N0.*3)),BER,'-b.','MarkerSize',10);
legend("uncoding");
title('BitErrorRate to Es/N0 in AWGN');
xlabel('Eb/N0(dB)');
ylabel('BitErrorRate(BER)');





function [symbol_in]=MapToSymbol(BinToDec,syNum)
    index=zeros(1,syNum);
    index(BinToDec==0)=15;
    index(BinToDec==1)=11;
    index(BinToDec==2)=3;
    index(BinToDec==3)=7;
    index(BinToDec==4)=16;
    index(BinToDec==5)=12;
    index(BinToDec==6)=4;
    index(BinToDec==7)=8;
    index(BinToDec==8)=14;
    index(BinToDec==9)=10;
    index(BinToDec==10)=2;
    index(BinToDec==11)=6;
    index(BinToDec==12)=13;
    index(BinToDec==13)=9;
    index(BinToDec==14)=1;
    index(BinToDec==15)=5;
    symbol_in=index;
end

function [v]=modulation(symbol_in,sx,sy,syNum)
    index=zeros(1,syNum);
    for i=1:16
        index(symbol_in==i)=sx(i);
    end
    index1=index';

    for i=1:16
        index(symbol_in==i)=sy(i);
    end
    index2=index';
    v=[index1 index2];
end


function [bit_out]=MapToBits(symbol_out,syNum)
    index=cell(1,syNum);
    index(symbol_out==1)={[1 1 1 0]};
    index(symbol_out==2)={[1 0 1 0]};
    index(symbol_out==3)={[0 0 1 0]};
    index(symbol_out==4)={[0 1 1 0]};
    index(symbol_out==5)={[1 1 1 1]};
    index(symbol_out==6)={[1 0 1 1]};
    index(symbol_out==7)={[0 0 1 1]};
    index(symbol_out==8)={[0 1 1 1]};
    index(symbol_out==9)={[1 1 0 1]};
    index(symbol_out==10)={[1 0 0 1]};
    index(symbol_out==11)={[0 0 0 1]};
    index(symbol_out==12)={[0 1 0 1]};
    index(symbol_out==13)={[1 1 0 0]};
    index(symbol_out==14)={[1 0 0 0]};
    index(symbol_out==15)={[0 0 0 0]};
    index(symbol_out==16)={[0 1 0 0]};
    bit_out=(reshape(cell2mat(index),[4,syNum]))'; 
end



