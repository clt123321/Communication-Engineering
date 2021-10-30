Es=2.5*d*d;%enger of signal power(d^2)
sigma=0.1:0.01:1.2;
N0=2*(sigma.^2);% En=noise RMS value^2=N0/2


x=10.*log10(Es./N0);
figure(1);
semilogy(x,SER,'-r.','MarkerSize',7);
hold on;
semilogy(x,uncoded_SER_theory,'-y.','MarkerSize',5);
hold on;
semilogy(x,SER_revised,'-b.','MarkerSize',7);
hold on;
semilogy(x,coded_SER_theory,'-g.','MarkerSize',5);
line([0,25],[10^-4,10^-4],'linestyle','--');
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
line([0,25],[10^-4,10^-4],'linestyle','--');
legend("uncoded BER","uncoded BER theory","coded BER","coded BER theory");
title('BitErrorRate to Es/N0 in AWGN');
xlabel('Eb/N0(dB)');
ylabel('BitErrorRate(BER)');

