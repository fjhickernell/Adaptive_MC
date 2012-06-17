%Time versus dimension
close all
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)

figure;
loglog(test.randch.dim,res.iidtime,'b.','markersize',20)
xlabel('{\it d}')
ylabel('Time (seconds)')
axis([1 100 0.001 100])
print -depsc iidtimedim.eps

figure;
loglog(test.randch.dim,res.iidheavytime,'b.','markersize',20)
xlabel('{\it d}')
ylabel('Time (seconds)')
axis([1 100 0.001 100])
print -depsc iidheavytimedim.eps

figure;
loglog(test.randch.dim,res.Soboltime,'b.','markersize',20)
xlabel('{\it d}')
ylabel('Time (seconds)')
axis([1 100 0.001 100])
print -depsc Soboltimedim.eps

figure;
loglog(test.randch.dim,res.Sobolheavytime,'b.','markersize',20)
xlabel('{\it d}')
ylabel('Time (seconds)')
axis([1 100 0.001 100])
print -depsc Sobolheavytimedim.eps

