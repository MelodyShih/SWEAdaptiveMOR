close all;
m = 500;

load('data/res1.mat');
load('data/restime1.mat')
t1 = time;
semilogy(ressqr(1:m),'--k');
hold on;

load('data/res2.mat');
load('data/restime2.mat')
t2 = time;
semilogy(ressqr(1:m),'--.k');

load('data/res5.mat');
load('data/restime5.mat')
t5 = time;
semilogy(ressqr(1:m),'-.k');

load('data/res10.mat');
load('data/restime10.mat')
t10 = time;
semilogy(ressqr(1:m),'-k');

legend({['$t =$',num2str(t1)],['$t =$',num2str(t2)], ['$t =$',num2str(t5)], ['$t =$',num2str(t10)]}, ...
       'Interpreter', 'latex','FontSize',15);
% title('Squared  $r_k = q_k - U_k(P_k^TU_k)^{-1}P_k^Tq_k$', 'Interpreter', 'latex','FontSize',15);
ylabel('residual square','Interpreter', 'latex','FontSize',15)
xlabel('index','Interpreter', 'latex','FontSize',15)
grid on;