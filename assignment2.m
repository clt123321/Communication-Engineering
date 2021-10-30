clear all;
clc
%
%17206018 ChengLitao
%tasks:
%1.plot show the simulated and theoretical bit error rate (SER) versusE b /N 0 in dB.
%2.plot show the simulated and theoretical bit error rate (BER) versusE b /N 0 in dB.
%3.find cross point which the Hamming codeoffers improved performance over an uncoded system(BER)
%4.estimate the value of E s /N 0 above which the system BER lies below 10^4
%**************************************************************************************
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
%**************************************************************************************


%**************************************************************************************
%simulation 16-QAM with hamming code
groupNum=150000;
symbolNum=4*groupNum;
d=1;
Es=2.5*d*d;%enger of signal power(d^2)

sx=[-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5];
sy=[1.5,1.5,1.5,1.5,0.5,0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5,-1.5,-1.5,-1.5,-1.5];
SER=[];
BER=[];
BER_revised=[];
SER_revised=[];
for sigma=0.1:0.01:1.2
    sigma
    %1.CreatbitSignal
    data_bits_in=rand(symbolNum,4)<0.5;
    %2.Hamming encode
    code_bits_in=mod(data_bits_in*G,2);
    bit_in=(reshape(code_bits_in,[4,7*groupNum]))';


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

    %6.Hamming decode
    bits_out=MapToBits(symbol_out,7*groupNum);
    code_bits_out=(reshape(bits_out',[symbolNum,7]));
    syndrome=mod(code_bits_out*H',2);
    errorVector= syndromToerror(syndrome,4*groupNum);
    revised_bits=mod(code_bits_out+errorVector,2);
    
    %7calculate SER & BER
    a1= (symbol_in'==symbol_out);
    SER=[SER,1-(sum(a1(:)==1)/(7*groupNum))];

    a3=bit_in==bits_out;
    BER=[BER,1-(sum(a3(:)==1)/(4*7*groupNum))];
    
    a4= (code_bits_in(:,4:7)==revised_bits(:,4:7));
    BER_revised=[BER_revised,1-(sum(a4(:)==1)/(4*4*groupNum))];
    
    a5= a4(:,1)&a4(:,2)&a4(:,3)&a4(:,4);
    SER_revised=[SER_revised,1-(sum(a5(:)==1)/(4*groupNum))];   
end
%**************************************************************************************
%
%Derivation of the theoretical value
%
%**************************************************************************************
sigma=0.1:0.01:1.2;
N0=2*(sigma.^2);% En=noise RMS value^2=N0/2 

gain =(2^4-1)*3*d^2/(6*Es*sqrt(2));
%P4=(2*(4-1)/4).*qfunc(d./sqrt(2.*N0));                       % no gain, 4-ASK SER_theory
P4=(2*(4-1)/4).*qfunc(sqrt(gain)*d./sqrt(2.*N0));   % coded_SER_theory = (2*(M-1)/M).*Q(sqrt(coding_gain.*d^2./2.*N0)), 4-ASK SER_theory
P4_bit =P4./2;                                       % coded_BER_theory = Pe/log2(4), 4-ASK BER_theory
uncoded_SER_theory=3*qfunc(d./(sqrt(2*N0)))-2.25.*(qfunc(d./(sqrt(2.*N0)))).^2;          % calculate uncoded 16-QAM SER_theory
coded_SER_theory = 1-(1-P4).^2;                                                          % calculate coded 16-QAM SER_theory


uncoded_BER_theory=3/8.*erfc(sqrt(0.4*0.625./N0))+0.25.*erfc(3.*sqrt(0.4*0.625./N0));    % calculate uncoded 16-QAM BER_theory
coded_BER_theory =1-(1-P4_bit).^2;                                                       % calculate coded 16-QAM BER_theory
%coded_BER_theory=3/8.*erfc(sqrt(0.4*0.625*gain./N0))+0.25.*erfc(3.*sqrt(0.4*0.625*gain./N0));  
%**************************************************************************************


%**************************************************************************************
%plot SER to Es/N0
x=10.*log10(Es./N0);
figure(1);
semilogy(x,SER,'-r.','MarkerSize',7);
hold on;
semilogy(x,uncoded_SER_theory,'-y.','MarkerSize',5);
hold on;
semilogy(x,SER_revised,'-b.','MarkerSize',7);
hold on;
semilogy(x,coded_SER_theory,'-g.','MarkerSize',5);
line([-5,25],[10^-4,10^-4],'linestyle','--');
legend("uncoded SER","uncoded SER theory","coded SER","coded SER theory");
title('SymbolErrorRate to Es/N0 in AWGN');
xlabel('Es/N0(dB)');
ylabel('SymbolErrorRate(SER)');
hold off;

%plot BER to Eb/N0
figure(2);
semilogy(10.*log10(Es./(N0.*4)),BER,'-r.','MarkerSize',7);
hold on;
semilogy(10.*log10(Es./(N0.*4)),uncoded_BER_theory,'-y.','MarkerSize',5);
hold on;
semilogy(10.*log10(Es./(N0.*4)),BER_revised,'-b.','MarkerSize',7);
hold on;
semilogy(10.*log10(Es./(N0.*4)),coded_BER_theory,'-g.','MarkerSize',5);
line([-10,25],[10^-4,10^-4],'linestyle','--');
legend("uncoded BER","uncoded BER theory","coded BER","coded BER theory");
title('BitErrorRate to Es/N0 in AWGN');
xlabel('Eb/N0(dB)');
ylabel('BitErrorRate(BER)');
%**************************************************************************************

%**************************************************************************************

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

function [errorVector]= syndromToerror(syndrome,syNum)
    BinToDec=syndrome(:,1)*4+syndrome(:,2)*2+syndrome(:,3);
    index=cell(1,syNum);
    index(BinToDec==0)={[0 0 0 0 0 0 0]};
    index(BinToDec==1)={[0 0 1 0 0 0 0]};
    index(BinToDec==2)={[0 1 0 0 0 0 0]};
    index(BinToDec==3)={[0 0 0 0 1 0 0]};
    index(BinToDec==4)={[1 0 0 0 0 0 0]};
    index(BinToDec==5)={[0 0 0 0 0 0 1]};
    index(BinToDec==6)={[0 0 0 1 0 0 0]};
    index(BinToDec==7)={[0 0 0 0 0 1 0]};
    errorVector=(reshape(cell2mat(index),[7,syNum]))';
end