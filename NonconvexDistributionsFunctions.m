%Picture of distributions and functions
%  whose convex combinations have higher kurtoses
clear all, close all
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)

%% Distributions
figure;
subplot(3,1,1)
x1=[-1 1;-1 1];
f1=[0 0;1 1]/2;
plot(x1,f1,'b-','linewidth',2)
axis([-3 3 0 1])
%xlabel('{\it y}')
ylabel('{\it f}_1')
text(-2.8,0.8,'\kappa=1')

subplot(3,1,2)
x2=[-1 1;-1 1]*2;
f2=[0 0;1 1]/2;
plot(x2,f2,'b-','linewidth',2)
axis([-3 3 0 1])
ylabel('{\it f}_2')
text(-2.8,0.8,'\kappa=1')

subplot(3,1,3)
x3=[-2 -1 1 2;-2 -1 1 2];
f3=[0 0 0 0;1 1 1 1]/4;
plot(x3,f3,'b-','linewidth',2)
axis([-3 3 0 1])
ylabel('{\it f}_3')
text(-2.8,0.8,'\kappa=34/25')
print -depsc ConvexSumDist.eps

%% Functions
figure;
subplot(3,1,1)
x1=[0 1/2 1/2 1];
f1=[-1 -1 1 1];
plot(x1,f1,'b-','linewidth',2)
axis([-0.1 1.1 -1.2 1.2])
%xlabel('{\it y}')
ylabel('{\it f}_1')
text(0,0.7,'\kappa=1')

subplot(3,1,2)
x2=[0 1/4 1/4 3/4 3/4 1];
f2=[1 1 -1 -1 1 1];
plot(x2,f2,'b-','linewidth',2)
axis([-0.1 1.1 -1.2 1.2])
ylabel('{\it f}_2')
text(0,0.7,'\kappa=1')

subplot(3,1,3)
x3=[0 1/4 1/4 1/2 1/2 3/4 3/4 1];
f3=[0 0 -1 -1 0 0 1 1];
plot(x3,f3,'b-','linewidth',2)
axis([-0.1 1.1 -1.2 1.2])
ylabel('{\it f}_3')
text(0,0.7,'\kappa=2')
print -depsc ConvexSumFun.eps


