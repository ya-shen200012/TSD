clear
ITER=10;
kappa = 1e2;
tola=[1e-6 1e-9 1e-12];
n = 1e4;
Maxiter = 1e5;
s1.iter=zeros(7,3);
%s1.time=zeros(7,3);
s1.time2=zeros(7,3);
s1.gnorm=zeros(7,3);
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
d(n/2+1:n)=1 +(kappa-1).*(0.2.*rand(n/2,1));
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
%A=diag(A);
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
s.iter=zeros(7,ITER+1);
s.gnorm=zeros(7,ITER+1);
%s.time=zeros(7,ITER+1);
s.time2=zeros(7,ITER+1);
for i=1:ITER
    disp(i);
x0 = -10 + 20*rand(n,1);
tol=tola(j);
maxabsgrad=max(abs(grad(x0)));
tol=tol* maxabsgrad;
disp('sd_21')
tic
[x1,iter,gnorm] = sd_2(x0, A,SD_2N1,  grad, tol, Maxiter);
s.iter(1,i)=iter;
%s.time(1,i)=out.cputime;
s.gnorm(1,i)=gnorm;
s.time2(1,i)=toc;
disp('sd_22')
tic
[x1,iter,gnorm] = sd_2(x0, A,SD_2N2,  grad, tol, Maxiter);
s.iter(2,i)=iter;
%s.time(2,i)=out.cputime;
s.gnorm(2,i)=gnorm;
s.time2(2,i)=toc;
disp('sd_23')
tic
[x1,iter,gnorm] = sd_2(x0, A,SD_2N3,  grad, tol, Maxiter);
s.iter(3,i)=iter;
%s.time(3,i)=out.cputime;
s.gnorm(3,i)=gnorm;
s.time2(3,i)=toc;
disp('sd_24')
tic
[x1,iter,gnorm] = sd_2(x0, A,SD_2N4,  grad, tol, Maxiter);
s.iter(4,i)=iter;
%s.time(4,i)=out.cputime;
s.gnorm(4,i)=gnorm;
s.time2(4,i)=toc;
disp('sd_25')
tic
[x1,iter,gnorm] = sd_2(x0, A,SD_2N5,  grad, tol, Maxiter);
s.iter(5,i)=iter;
%s.time(5,i)=out.cputime;
s.gnorm(5,i)=gnorm;
s.time2(5,i)=toc;
disp('sd_26')
tic
[x1,iter,gnorm] = sd_2(x0, A,SD_2N6,  grad, tol, Maxiter);
s.iter(6,i)=iter;
%s.time(6,i)=out.cputime;
s.gnorm(6,i)=gnorm;
s.time2(6,i)=toc;
disp('sd_27')
tic
[x1,iter,gnorm] = sd_2(x0, A,SD_2N7,  grad, tol,Maxiter);
s.iter(7,i)=iter;
%s.time(7,i)=out.cputime;
s.gnorm(7,i)=gnorm;
s.time2(7,i)=toc;
% opts.PrintLevel = 0;
% opts.tol=tol;
% opts.tau=0.9;
% opts.gama=1.05;
% disp('bbq')
% tic
% [out, hist] = BBQ_unc(x0, fobj, grad, opts);
% s.iter(19,i)=out.iter;
% %s(1).time(19,i)=out.cputime;
% %s.gnorm(19,i)=out.gnorm;
% s.time2(19,i)=toc;
end
for i=1:7
s.iter(i,ITER+1)=mean(s.iter(i,1:ITER));
%s.time(i,ITER+1)=mean(s.time(i,1:ITER));
s.time2(i,ITER+1)=mean(s.time2(i,1:ITER));
s.gnorm(i,ITER+1)=mean(s.gnorm(i,1:ITER));
end
s1.iter(1,j)=s.iter(1,ITER+1);
s1.iter(2,j)=s.iter(2,ITER+1);
s1.iter(3,j)=s.iter(3,ITER+1);
s1.iter(4,j)=s.iter(4,ITER+1);
s1.iter(5,j)=s.iter(5,ITER+1);
s1.iter(6,j)=s.iter(6,ITER+1);
s1.iter(7,j)=s.iter(7,ITER+1);
%s1.iter(8,j)=s.iter(18,ITER+1);
s1.time2(1,j)=s.time2(1,ITER+1);
s1.time2(2,j)=s.time2(2,ITER+1);
s1.time2(3,j)=s.time2(3,ITER+1);
s1.time2(4,j)=s.time2(4,ITER+1);
s1.time2(5,j)=s.time2(5,ITER+1);
s1.time2(6,j)=s.time2(6,ITER+1);
s1.time2(7,j)=s.time2(7,ITER+1);
s1.gnorm(1,j)=s.gnorm(1,ITER+1);
s1.gnorm(2,j)=s.gnorm(2,ITER+1);
s1.gnorm(3,j)=s.gnorm(3,ITER+1);
s1.gnorm(4,j)=s.gnorm(4,ITER+1);
s1.gnorm(5,j)=s.gnorm(5,ITER+1);
s1.gnorm(6,j)=s.gnorm(6,ITER+1);
s1.gnorm(7,j)=s.gnorm(7,ITER+1);
end
