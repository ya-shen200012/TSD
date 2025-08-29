clear
ITER=10;
kappa = 1e2;
tola=[1e-6 1e-9 1e-12];
n = 1e4;
Maxiter = 5e4;
s1.iter=zeros(1,3);
s1.gnorm=zeros(1,3);
s1.time2=zeros(1,3);
s2.iter=zeros(1,3);
s2.gnorm=zeros(1,3);
s2.time2=zeros(1,3);
SD_2N1=10;
SD_2N2=50;
SD_2N3=100;
SD_2N4=200;
SD_2N5=300;
SD_2N6=500;
SD_2N7=1000;
rng(20230629)
% % set1
% d = 1 + (kappa-1).*rand(n,1);
% d(1) = 1; d(end) = kappa;
% set2
d=zeros(n,1);
d(1:n/2) = 1 +(kappa-1).*(0.8+0.2.*rand(n/2,1));
d((n/2+1):n)=1 +(kappa-1).*(0.2.*rand(n/2,1));
d(1) = 1; d(end) = kappa;
% % set3
% d=zeros(n,1);
% d(1:n/5) = 1 + (100-1).*rand(n/5,1);
% d(n/5+1:n)=kappa/2+(kappa/2).*rand(n/5*4,1);
% d(1) = 1; d(end) = kappa;
% % set4
% d=zeros(n,1);
% for i=1:n
% d(i)=kappa^((n-i)/(n-1));
% end
% % set5
% d=zeros(n,1);
% d(1:n*4/5) = 1 + (100-1).*rand(n*4/5,1);
% d(n*4/5+1:n)=kappa/2+(kappa/2).*rand(n/5,1);
% d(1) = 1; d(end) = kappa;
A = d(:);
% A=diag(A);
xs = -10 + 20*rand(n,1);

if size(A,2) == 1
    grad = @(x) A.*(x - xs);
    fobj = @(x) (x - xs)'*(A.*(x - xs));
    %fAx = @(x) A.*x;
else
    grad =   @(x) A*(x - xs);
    fobj = @(x) (x - xs)'*(A*(x - xs));
    %fAx = @(x) A*x;
end
for j=1:3
s.iter=zeros(54,ITER+1);
s.gnorm=zeros(54,ITER+1);
s.time2=zeros(54,ITER+1);
for i=1:ITER
    disp(i);
x0 = -10 + 20*rand(n,1);
tol=tola(j);
maxabsgrad=max(abs(grad(x0)));
tol=tol* maxabsgrad;
opts.PrintLevel = 1;
opts.tol=tol;
opts.gama=1;
opts.tau=0.1;
disp('bbq1')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(1,i)=out.iter;
%s(1).time(1,i)=out.cputime;
s.gnorm(1,i)=out.gnorm;
s.time2(1,i)=toc;
opts.tau=0.2;
disp('bbq2')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(2,i)=out.iter;
%s(1).time(2,i)=out.cputime;
s.gnorm(2,i)=out.gnorm;
s.time2(2,i)=toc;
opts.tau=0.3;
disp('bbq3')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(3,i)=out.iter;
%s(1).time(3,i)=out.cputime;
s.gnorm(3,i)=out.gnorm;
s.time2(3,i)=toc;
opts.tau=0.4;
disp('bbq4')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(4,i)=out.iter;
%s(1).time(4,i)=out.cputime;
s.gnorm(4,i)=out.gnorm;
s.time2(4,i)=toc;
opts.tau=0.5;
disp('bbq5')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(5,i)=out.iter;
%s(1).time(5,i)=out.cputime;
s.gnorm(5,i)=out.gnorm;
s.time2(5,i)=toc;
opts.tau=0.6;
disp('bbq6')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(6,i)=out.iter;
%s(1).time(6,i)=out.cputime;
s.gnorm(6,i)=out.gnorm;
s.time2(6,i)=toc;
opts.tau=0.7;
disp('bbq7')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(7,i)=out.iter;
%s(1).time(7,i)=out.cputime;
s.gnorm(7,i)=out.gnorm;
s.time2(7,i)=toc;
opts.tau=0.8;
disp('bbq8')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(8,i)=out.iter;
%s(1).time(8,i)=out.cputime;
s.gnorm(8,i)=out.gnorm;
s.time2(8,i)=toc;
opts.tau=0.9;
disp('bbq9')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(9,i)=out.iter;
%s(1).time(9,i)=out.cputime;
s.gnorm(9,i)=out.gnorm;
s.time2(9,i)=toc;

%222
opts.gama=1.02;
opts.tau=0.1;
disp('bbq10')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(10,i)=out.iter;
%s(1).time(10,i)=out.cputime;
s.gnorm(10,i)=out.gnorm;
s.time2(10,i)=toc;
opts.tau=0.2;
disp('bbq11')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(11,i)=out.iter;
%s(1).time(11,i)=out.cputime;
s.gnorm(11,i)=out.gnorm;
s.time2(11,i)=toc;
opts.tau=0.3;
disp('bbq12')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(12,i)=out.iter;
%s(1).time(12,i)=out.cputime;
s.gnorm(12,i)=out.gnorm;
s.time2(12,i)=toc;
opts.tau=0.4;
disp('bbq13')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(13,i)=out.iter;
%s(1).time(13,i)=out.cputime;
s.gnorm(13,i)=out.gnorm;
s.time2(13,i)=toc;
opts.tau=0.5;
disp('bbq14')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(14,i)=out.iter;
%s(1).time(14,i)=out.cputime;
s.gnorm(14,i)=out.gnorm;
s.time2(14,i)=toc;
opts.tau=0.6;
disp('bbq15')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(15,i)=out.iter;
%s(1).time(15,i)=out.cputime;
s.gnorm(15,i)=out.gnorm;
s.time2(15,i)=toc;
opts.tau=0.7;
disp('bbq16')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(16,i)=out.iter;
%s(1).time(16,i)=out.cputime;
s.gnorm(16,i)=out.gnorm;
s.time2(16,i)=toc;
opts.tau=0.8;
disp('bbq17')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(17,i)=out.iter;
%s(1).time(17,i)=out.cputime;
s.gnorm(17,i)=out.gnorm;
s.time2(17,i)=toc;
opts.tau=0.9;
disp('bbq18')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(18,i)=out.iter;
%s(1).time(18,i)=out.cputime;
s.gnorm(18,i)=out.gnorm;
s.time2(18,i)=toc;

%333
opts.gama=1.05;
opts.tau=0.1;
disp('bbq19')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(19,i)=out.iter;
%s(1).time(19,i)=out.cputime;
s.gnorm(19,i)=out.gnorm;
s.time2(19,i)=toc;
opts.tau=0.2;
disp('bbq20')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(20,i)=out.iter;
%s(1).time(20,i)=out.cputime;
s.gnorm(20,i)=out.gnorm;
s.time2(20,i)=toc;
opts.tau=0.3;
disp('bbq21')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(21,i)=out.iter;
%s(1).time(21,i)=out.cputime;
s.gnorm(21,i)=out.gnorm;
s.time2(21,i)=toc;
opts.tau=0.4;
disp('bbq22')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(22,i)=out.iter;
%s(1).time(22,i)=out.cputime;
s.gnorm(22,i)=out.gnorm;
s.time2(22,i)=toc;
opts.tau=0.5;
disp('bbq23')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(23,i)=out.iter;
%s(1).time(23,i)=out.cputime;
s.gnorm(23,i)=out.gnorm;
s.time2(23,i)=toc;
opts.tau=0.6;
disp('bbq24')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(24,i)=out.iter;
%s(1).time(24,i)=out.cputime;
s.gnorm(24,i)=out.gnorm;
s.time2(24,i)=toc;
opts.tau=0.7;
disp('bbq25')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(25,i)=out.iter;
%s(1).time(25,i)=out.cputime;
s.gnorm(25,i)=out.gnorm;
s.time2(25,i)=toc;
opts.tau=0.8;
disp('bbq26')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(26,i)=out.iter;
%s(1).time(26,i)=out.cputime;
s.gnorm(26,i)=out.gnorm;
s.time2(26,i)=toc;
opts.tau=0.9;
disp('bbq27')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(27,i)=out.iter;
%s(1).time(27,i)=out.cputime;
s.gnorm(27,i)=out.gnorm;
s.time2(27,i)=toc;

%444
opts.gama=1.1;
opts.tau=0.1;
disp('bbq28')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(28,i)=out.iter;
%s(1).time(28,i)=out.cputime;
s.gnorm(28,i)=out.gnorm;
s.time2(28,i)=toc;
opts.tau=0.2;
disp('bbq29')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(29,i)=out.iter;
%s(1).time(29,i)=out.cputime;
s.gnorm(29,i)=out.gnorm;
s.time2(29,i)=toc;
opts.tau=0.3;
disp('bbq30')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(30,i)=out.iter;
%s(1).time(30,i)=out.cputime;
s.gnorm(30,i)=out.gnorm;
s.time2(30,i)=toc;
opts.tau=0.4;
disp('bbq31')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(31,i)=out.iter;
%s(1).time(31,i)=out.cputime;
s.gnorm(31,i)=out.gnorm;
s.time2(31,i)=toc;
opts.tau=0.5;
disp('bbq32')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(32,i)=out.iter;
%s(1).time(32,i)=out.cputime;
s.gnorm(32,i)=out.gnorm;
s.time2(32,i)=toc;
opts.tau=0.6;
disp('bbq33')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(33,i)=out.iter;
%s(1).time(33,i)=out.cputime;
s.gnorm(33,i)=out.gnorm;
s.time2(33,i)=toc;
opts.tau=0.7;
disp('bbq34')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(34,i)=out.iter;
%s(1).time(34,i)=out.cputime;
s.gnorm(34,i)=out.gnorm;
s.time2(34,i)=toc;
opts.tau=0.8;
disp('bbq35')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(35,i)=out.iter;
%s(1).time(35,i)=out.cputime;
s.gnorm(35,i)=out.gnorm;
s.time2(35,i)=toc;
opts.tau=0.9;
disp('bbq36')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(36,i)=out.iter;
%s(1).time(36,i)=out.cputime;
s.gnorm(36,i)=out.gnorm;
s.time2(36,i)=toc;

%555
opts.gama=1.2;
opts.tau=0.1;
disp('bbq37')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(37,i)=out.iter;
%s(1).time(37,i)=out.cputime;
s.gnorm(37,i)=out.gnorm;
s.time2(37,i)=toc;
opts.tau=0.2;
disp('bbq38')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(38,i)=out.iter;
%s(1).time(38,i)=out.cputime;
s.gnorm(38,i)=out.gnorm;
s.time2(38,i)=toc;
opts.tau=0.3;
disp('bbq39')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(39,i)=out.iter;
%s(1).time(39,i)=out.cputime;
s.gnorm(39,i)=out.gnorm;
s.time2(39,i)=toc;
opts.tau=0.4;
disp('bbq40')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(40,i)=out.iter;
%s(1).time(40,i)=out.cputime;
s.gnorm(40,i)=out.gnorm;
s.time2(40,i)=toc;
opts.tau=0.5;
disp('bbq41')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(41,i)=out.iter;
%s(1).time(41,i)=out.cputime;
s.gnorm(41,i)=out.gnorm;
s.time2(41,i)=toc;
opts.tau=0.6;
disp('bbq42')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(42,i)=out.iter;
%s(1).time(42,i)=out.cputime;
s.gnorm(42,i)=out.gnorm;
s.time2(42,i)=toc;
opts.tau=0.7;
disp('bbq43')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(43,i)=out.iter;
%s(1).time(43,i)=out.cputime;
s.gnorm(43,i)=out.gnorm;
s.time2(43,i)=toc;
opts.tau=0.8;
disp('bbq44')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(44,i)=out.iter;
%s(1).time(44,i)=out.cputime;
s.gnorm(44,i)=out.gnorm;
s.time2(44,i)=toc;
opts.tau=0.9;
disp('bbq45')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(45,i)=out.iter;
%s(1).time(45,i)=out.cputime;
s.gnorm(45,i)=out.gnorm;
s.time2(45,i)=toc;
% 666 
opts.gama=1.3;
opts.tau=0.1;
disp('bbq46')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(46,i)=out.iter;
%s(1).time(46,i)=out.cputime;
s.gnorm(46,i)=out.gnorm;
s.time2(46,i)=toc;
opts.tau=0.2;
disp('bbq47')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(47,i)=out.iter;
%s(1).time(47,i)=out.cputime;
s.gnorm(47,i)=out.gnorm;
s.time2(47,i)=toc;
opts.tau=0.3;
disp('bbq48')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(48,i)=out.iter;
%s(1).time(48,i)=out.cputime;
s.gnorm(48,i)=out.gnorm;
s.time2(48,i)=toc;
opts.tau=0.4;
disp('bbq49')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(49,i)=out.iter;
%s(1).time(49,i)=out.cputime;
s.gnorm(49,i)=out.gnorm;
s.time2(49,i)=toc;
opts.tau=0.5;
disp('bbq50')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(50,i)=out.iter;
%s(1).time(50,i)=out.cputime;
s.gnorm(50,i)=out.gnorm;
s.time2(50,i)=toc;
opts.tau=0.6;
disp('bbq51')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(51,i)=out.iter;
%s(1).time(51,i)=out.cputime;
s.gnorm(51,i)=out.gnorm;
s.time2(51,i)=toc;
opts.tau=0.7;
disp('bbq52')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(52,i)=out.iter;
%s(1).time(52,i)=out.cputime;
s.gnorm(52,i)=out.gnorm;
s.time2(52,i)=toc;
opts.tau=0.8;
disp('bbq53')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(53,i)=out.iter;
%s(1).time(53,i)=out.cputime;
s.gnorm(53,i)=out.gnorm;
s.time2(53,i)=toc;
opts.tau=0.9;
disp('bbq54')
tic
[out, ~] = BBQ_unc(x0, fobj, grad, opts);
s.iter(54,i)=out.iter;
%s(1).time(54,i)=out.cputime;
s.gnorm(54,i)=out.gnorm;
s.time2(54,i)=toc;
end
for i=1:54
s.iter(i,ITER+1)=mean(s.iter(i,1:ITER));
%s.time(i,ITER+1)=mean(s.time(i,1:ITER));
s.gnorm(i,ITER+1)=mean(s.gnorm(i,1:ITER));
s.time2(i,ITER+1)=mean(s.time2(i,1:ITER));
end
% [s1.iter(1,j), s2.iter(1,j)]=min(s.iter(1:54,ITER+1));
% [s1.time2(1,j), s2.time2(1,j)]=min(s.time2(1:54,ITER+1));
[s1.iter(1,j), s2.iter(1,j)]=min(s.iter(1:54,ITER+1));
[s1.gnorm(1,j), s2.gnorm(1,j)]=max(s.gnorm(1:54,ITER+1));
[s1.time2(1,j), s2.time2(1,j)]=min(s.time2(1:54,ITER+1));
end
