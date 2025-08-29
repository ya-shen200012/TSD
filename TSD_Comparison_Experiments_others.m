clear
ITER=10;
kappa = 1e2;
tola=[1e-6 1e-9 1e-12];
n = 1e4;
Maxiter = 5e4;
s1.iter=zeros(3,3);
%s1.time=zeros(3,3);
s1.time2=zeros(3,3);
s1.gnorm=zeros(3,3);
s2.iter=zeros(1,3);
%s2.time=zeros(1,3);
s2.time2=zeros(1,3);
s2.gnorm=zeros(1,3);
abbmin2tao1=0.1;
abbmin2tao2=0.2;
abbmin2tao3=0.3;
abbmin2tao4=0.4;
abbmin2tao5=0.5;
abbmin2tao6=0.6;
abbmin2tao7=0.7;
abbmin2tao8=0.8;
abbmin2tao9=0.9;
rng(20230629)
% set1
d = 1 + (kappa-1).*rand(n,1);
d(1) = 1; d(end) = kappa;
% % set2
% d=zeros(n,1);
% d(1:n/2) = 1 +(kappa-1).*(0.8+0.2.*rand(n/2,1));
% d(n/2+1:n)=1 +(kappa-1).*(0.2.*rand(n/2,1));
% d(1) = 1; d(end) = kappa;
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
s.iter=zeros(22,ITER+1);
%s.time=zeros(22,ITER+1);
s.time2=zeros(22,ITER+1);
s.gnorm=zeros(22,ITER+1);
for i=1:ITER
    disp(i);
x0 = -10 + 20*rand(n,1);
tol=tola(j);
maxabsgrad=max(abs(grad(x0)));
tol=tol* maxabsgrad;
disp('BB1')
tic
[~,iter,gnorm] = bb1(x0, A,  grad, tol);
s.iter(1,i)=iter;
%s.time(1,i)=out.cputime;
s.gnorm(1,i)=gnorm;
s.time2(1,i)=toc;
disp('dy')
tic
[~,iter,gnorm] = dy(x0, A,  grad, tol);
s.iter(2,i)=iter;
%s.time(2,i)=out.cputime;
s.gnorm(2,i)=gnorm;
s.time2(2,i)=toc;
disp('abbmin21')
tic
[~,iter,gnorm] = abbmin2(x0, A, grad,tol,abbmin2tao1);
s.iter(3,i)=iter;
%s.time(3,i)=out.cputime;
s.gnorm(3,i)=gnorm;
s.time2(3,i)=toc;
disp('abbmin22')
tic
[~,iter,gnorm] = abbmin2(x0, A, grad,tol,abbmin2tao2);
s.iter(4,i)=iter;
%s.time(4,i)=out.cputime;
s.gnorm(4,i)=gnorm;
s.time2(4,i)=toc;
disp('abbmin23')
tic
[~,iter,gnorm] = abbmin2(x0, A, grad,tol,abbmin2tao3);
s.iter(5,i)=iter;
%s.time(5,i)=out.cputime;
s.gnorm(5,i)=gnorm;
s.time2(5,i)=toc;
disp('abbmin24')
tic
[~,iter,gnorm] = abbmin2(x0, A, grad,tol,abbmin2tao4);
s.iter(6,i)=iter;
%s.time(6,i)=out.cputime;
s.gnorm(6,i)=gnorm;
s.time2(6,i)=toc;
disp('abbmin25')
tic
[~,iter,gnorm] = abbmin2(x0, A, grad,tol,abbmin2tao5);
s.iter(7,i)=iter;
%s.time(7,i)=out.cputime;
s.gnorm(7,i)=gnorm;
s.time2(7,i)=toc;
disp('abbmin26')
tic
[~,iter,gnorm] = abbmin2(x0, A, grad,tol,abbmin2tao6);
s.iter(8,i)=iter;
%s.time(8,i)=out.cputime;
s.gnorm(8,i)=gnorm;
s.time2(8,i)=toc;
disp('abbmin27')
tic
[~,iter,gnorm] = abbmin2(x0, A, grad,tol,abbmin2tao7);
s.iter(9,i)=iter;
%s.time(9,i)=out.cputime;
s.gnorm(9,i)=gnorm;
s.time2(9,i)=toc;
disp('abbmin28')
tic
[~,iter,gnorm] = abbmin2(x0, A, grad,tol,abbmin2tao8);
s.iter(10,i)=iter;
%s.time(10,i)=out.cputime;
s.gnorm(10,i)=gnorm;
s.time2(10,i)=toc;
disp('abbmin29')
tic
[~,iter,gnorm] = abbmin2(x0, A, grad,tol,abbmin2tao9);
s.iter(11,i)=iter;
%s.time（11,i)=out.cputime;
s.gnorm(11,i)=gnorm;
s.time2(11,i)=toc;
end
for i=1:11
s.iter(i,ITER+1)=mean(s.iter(i,1:ITER));
%s.time(i,ITER+1)=mean(s.time(i,1:ITER));
s.time2(i,ITER+1)=mean(s.time2(i,1:ITER));
s.gnorm(i,ITER+1)=mean(s.gnorm(i,1:ITER));
end
s1.iter(1,j)=s.iter(1,ITER+1);
s1.iter(2,j)=s.iter(2,ITER+1);
[s1.iter(3,j), s2.iter(1,j)]=min(s.iter(3:11,ITER+1));
s1.time2(1,j)=s.time2(1,ITER+1);
s1.time2(2,j)=s.time2(2,ITER+1);
[s1.time2(3,j), s2.time2(1,j)]=min(s.time2(3:11,ITER+1));
s1.gnorm(1,j)=s.gnorm(1,ITER+1);
s1.gnorm(2,j)=s.gnorm(2,ITER+1);
[s1.gnorm(3,j), s2.gnorm(1,j)]=max(s.gnorm(3:11,ITER+1));
end
