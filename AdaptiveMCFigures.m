% Plot the output of the adaptive MC Example
set(0,'defaultaxesfontsize',18,'defaulttextfontsize',18)
clear all, close all

load AdaptiveMCOut.mat

whshow=1:6;
nshow=size(whshow,2);
legtext=cell(nshow);
for ii=1:nshow
    legtext{ii}=num2str(pvec(ii),'%3.0g');
end
probvec=((1:nrep)'-1/2)/nrep;
figure;
h=plot(sort(errBC(:,whshow)/epsilon),probvec,'linewidth',2);
legend(h,legtext{whshow},'location','northwest')
set(gca,'XScale','log')
axis([1e-4 1e2 0 1])
xlabel('Normalized Error')
ylabel('Probability')
print -depsc NormalErrFig.eps

kurtvec=1./(pvec.*(1-pvec))-3;

fprintf(1,'%8.4g & ',pvec(whshow)); fprintf(1,'\\\\ \r');
fprintf(1,'%8.2f & ',kurtvec(whshow)); fprintf(1,'\\\\ \r');
fprintf(1,'%8.2f\\%% & ',100*probsucc(whshow)); fprintf(1,'\\\\ \r');

% figure;
% h=plot(sort(sig./sig0up),probvec,'linewidth',2);
% set(gca,'XScale','log')
% axis([1e-2 1e0 0 1])
% legend(h,{'0.001','0.002','0.005','0.01','0.02','0.05','0.1','0.2','0.5'},'location','northwest')
