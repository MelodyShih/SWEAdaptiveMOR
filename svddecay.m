close all;
m = 4;

load('data/S0_30.mat');
semilogy(s,'-ko');
xlim([1,30])
hold on;

load('data/S1_30.mat');
semilogy(s,'-ks');

load('data/S2_30.mat');
semilogy(s,'-.ko');

load('data/S5_30.mat');
semilogy(s,'-k*');

load('data/S10_30.mat');
semilogy(s,'-kd');

load('data/Sall.mat');
semilogy(s(1:30),'-k.','MarkerSize',10);

legend({'from $t = 0$','from $t = 1$', 'from $t = 2$', 'from $t = 5$', 'from $t = 10$', 'global'}, ...
       'Interpreter', 'latex','FontSize',12);
title('Time steps $w= 30$', 'Interpreter', 'latex','FontSize',15);
ylabel('singular values','Interpreter', 'latex','FontSize',15)
xlabel('index','Interpreter', 'latex','FontSize',15)
grid on;