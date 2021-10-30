clear all;
clc
%
%17206018 ChengLitao
%
%
% plot (SER) versus E s /N 0 curve   SER on a log scale and E s /N 0 in dB.
% in the same graph, plot the theoretical SER
% On a new graph, plot the bit error rate (BER) versus E b /N 0 
% estimate the value of E b /N 0 above which the system BER lies below 10^4

symbolNum=100
d=1;
Es=3*d*d/2;%enger of signal power(d^2)
sx=[-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5];
sy=[1.5,1.5,0.5,0.5,-0.5,-0.5,-1.5,-1.5];
SER=[];
BER=[];

for sigma=0.1:0.01:0.5
    %1CreatbitSignal
    bit_in=rand(symbolNum,3)<0.5;

    %2MapToSymbol
    BinToDec=bit_in(:,1)*4+bit_in(:,2)*2+bit_in(:,3);
    index=zeros(1,symbolNum);
    index(BinToDec==0)=2;
    index(BinToDec==1)=4;
    index(BinToDec==2)=8;
    index(BinToDec==3)=6;
    index(BinToDec==4)=1;
    index(BinToDec==5)=3;
    index(BinToDec==6)=7;
    index(BinToDec==7)=5;
    symbol_in=index;


    %3modulate
    index(symbol_in==1)=sx(1);
    index(symbol_in==2)=sx(2);
    index(symbol_in==3)=sx(3);
    index(symbol_in==4)=sx(4);
    index(symbol_in==5)=sx(5);
    index(symbol_in==6)=sx(6);
    index(symbol_in==7)=sx(7);
    index(symbol_in==8)=sx(8);
    index1=index';

    index(symbol_in==1)=sy(1);
    index(symbol_in==2)=sy(2);
    index(symbol_in==3)=sy(3);
    index(symbol_in==4)=sy(4);
    index(symbol_in==5)=sy(5);
    index(symbol_in==6)=sy(6);
    index(symbol_in==7)=sy(7);
    index(symbol_in==8)=sy(8);
    index2=index';

    v=[index1 index2];

    %4 add noise
    noise=sigma.*randn(symbolNum,2);
    v=v+noise;

    %5 demodulate£¨optimum£©
    s=[sx' sy'];
    ds=zeros(symbolNum,8);
    for i=1:8
        ds(:,i)=(v(:,1)-sx(i)).^2+(v(:,2)-sy(i)).^2;
    end
    [min_d,symbol_out]=min(ds,[],2);
    symbol_out;

    %6 MapToBits
    index(symbol_out==1)=1;
    index(symbol_out==2)=0;
    index(symbol_out==3)=1;
    index(symbol_out==4)=0;
    index(symbol_out==5)=1;
    index(symbol_out==6)=0;
    index(symbol_out==7)=1;
    index(symbol_out==8)=0;
    b1=index';
    index(symbol_out==1)=0;
    index(symbol_out==2)=0;
    index(symbol_out==3)=0;
    index(symbol_out==4)=0;
    index(symbol_out==5)=1;
    index(symbol_out==6)=1;
    index(symbol_out==7)=1;
    index(symbol_out==8)=1;
    b2=index';
    index(symbol_out==1)=0;
    index(symbol_out==2)=0;
    index(symbol_out==3)=1;
    index(symbol_out==4)=1;
    index(symbol_out==5)=1;
    index(symbol_out==6)=1;
    index(symbol_out==7)=0;
    index(symbol_out==8)=0;
    b3=index';

    bits_out=[b1 b2 b3];

    %7calculate SER & BER
    a= (symbol_in'==symbol_out);
    SER=[SER,1-(sum(a(:)==1)/symbolNum)]


    b=bit_in==bits_out;
    BER=[BER,1-(sum(b(:)==1)/(3*symbolNum))]

end


sigma=0.1:0.01:0.5;
N0=2*(sigma.^2);% En=noise RMS value^2=N0/2
theory_value=2.5*qfunc(d./(sqrt(2*N0)))-1.5.*(qfunc(d./(sqrt(2.*N0)))).^2;

%plot SER to Es/N0
x=10.*log10(Es./N0);
figure(1);
semilogy(x,SER,'-b.','MarkerSize',10);
hold on;
semilogy(x,theory_value,'-r.','MarkerSize',5);
title('SymbolErrorRate to Es/N0 in AWGN');
xlabel('Es/N0(dB)');
ylabel('SymbolErrorRate(SER)');
hold off;


%plot BER to Eb/N0
figure(2);
semilogy(10.*log10(Es./(N0.*3)),BER,'-b.','MarkerSize',10);
title('BitErrorRate to Es/N0 in AWGN');
xlabel('Eb/N0(dB)');
ylabel('BitErrorRate(BER)');