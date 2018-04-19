%script to process data and make plots for DK paper

fsize=22;
pstn=[.2 .2 .4 .6];

% first example:  Compare DK with NSSV estimator
% u=exponential peak (w/eps^2=10^(-6)) + 0.01 sin(100 pi x) sin(100 pi y)
% Files  compnssv_estdk.mat, compnssv_estnssv.mat
load compnssv_estdk_new.mat

N_dk=errinfo(:,1);
err_dk=errinfo(:,2);
est_dk=errinfo(:,3)+errinfo(:,4)+errinfo(:,5);

load compnssv_estnssv_new.mat
N_nssv=errinfo(:,1);
err_nssv=errinfo(:,2);
est_nssv=errinfo(:,6)+errinfo(:,7)+errinfo(:,8);

hFig=figure(1)
set(hFig, 'units','normalized','Position', pstn)

loglog(N_nssv, est_nssv,'b-d','LineWidth',1.5,'MarkerSize',8)
hold on
loglog(N_nssv, err_nssv,'g-*','LineWidth',1.5,'MarkerSize',8)
loglog(N_dk, 150*log(N_dk)./N_dk,'k--','LineWidth',1.5, 'MarkerSize',8)
loglog(N_dk, est_dk,'r-o','LineWidth',1.5,'MarkerSize',8)
loglog(N_dk, err_dk,'c-+','LineWidth',1.5,'MarkerSize',8)
axis([10 10^7 10^(-4) 10^5])

set(gca,'FontSize',fsize)
set(gcf, 'Color', 'w');

 [hleg1, hobj1]=legend('$est_{NSSV}$' ,'$error_{NSSV}$',...
     '$\log(DOF)/DOF$','$est_{DK}$','$err_{DK}$');
 h=legend;
 set(h, 'interpreter', 'latex','Fontsize',fsize);
 set(hleg1,'units','normalized','position',[0 0 .38 .1])
 
 
 xlabel('DOF','interpreter', 'latex','Fontsize',fsize);
 ylabel('error','interpreter', 'latex','Fontsize',fsize);


%load testcf_1.mat
load newer_testcf1e-6_cf1e-0.mat
N_1=errinfo(:,1);
err_1=errinfo(:,2);
est_1=errinfo(:,3)+errinfo(:,4)+errinfo(:,5);
eff_inf_1=max(err_1./est_1)


%load testcf_.01.mat
load newer_testcf1e-6_cf1e-2.mat
N_01=errinfo(:,1);
lg=length(N_01);
N_01=N_01(1:lg-1);
err_01=errinfo(1:lg-1,2);
est_01=errinfo(1:lg-1,3)+errinfo(1:lg-1,4)+errinfo(1:lg-1,5);
eff_inf_01=max(err_01./est_01)




%load testcf_.0001.mat
load newer_testcf1e-6_cf1e-4.mat
N_0001=errinfo(:,1);
err_0001=errinfo(:,2);
est_0001=errinfo(:,3)+errinfo(:,4)+errinfo(:,5);
eff_inf_0001=max(err_0001./est_0001)



%load testcf_.000001.mat
load newer_testcf1e-6_cf1e-6.mat
N_000001=errinfo(:,1);
err_000001=errinfo(:,2);
est_000001=errinfo(:,3)+errinfo(:,4)+errinfo(:,5);

eff_inf_000001=max(err_000001./est_000001)

load testcf_.00000001.mat
N_00000001=errinfo(:,1);
err_00000001=errinfo(:,2);
est_00000001=errinfo(:,3)+errinfo(:,4)+errinfo(:,5);

eff_inf_00000001=max(err_00000001./est_00000001)


hFig=figure(2)
set(hFig, 'units','normalized','Position', pstn)

loglog(N_1, err_1,'k-d','LineWidth',1.5,'MarkerSize',8)
hold on
loglog(N_01, err_01,'k-*','LineWidth',1.5,'MarkerSize',8)
loglog(N_1, 150*log(N_1)./N_1,'k--','LineWidth',1.5, 'MarkerSize',8)
loglog(N_0001, err_0001,'k-o','LineWidth',1.5,'MarkerSize',8)
loglog(N_000001, err_000001,'k-+','LineWidth',1.5,'MarkerSize',8)
%loglog(N_1, est_1, 'k-x','LineWidth',1.5,'MarkerSize',8)


%loglog(N_00000001,err_00000001,'k-v','LineWidth',1.5,'MarkerSize',8)
axis([10 10^7 10^(-4) 1])

set(gca,'FontSize',fsize)
set(gcf, 'Color', 'w');

 [hleg1, hobj1]=legend('$C_f=1$','$C_f=10^{-2}$',...
     '$\log(DOF)/DOF$','$C_f=10^{-4}$','$C_f=10^{-6}$' );
 h=legend;
 set(h, 'interpreter', 'latex','Fontsize',fsize);
 set(hleg1,'units','normalized','position',[0 0 .39 .1])
 set(gca,'XTick',[10 10^3 10^5 10^7])
 
 xlabel('DOF','interpreter', 'latex','Fontsize',fsize);
 ylabel('error','interpreter', 'latex','Fontsize',fsize);

 
 
hFig=figure(100)
set(hFig, 'units','normalized','Position', pstn)
loglog(N_000001, (err_000001./est_000001).^(-1),'k-+','LineWidth',1.5,'MarkerSize',8)
hold on
loglog(N_0001, (err_0001./est_0001).^(-1),'k-o','LineWidth',1.5,'MarkerSize',8)
loglog(N_01, (err_01./est_01).^(-1),'k-*','LineWidth',1.5,'MarkerSize',8)
loglog(N_1, (err_1./est_1).^(-1),'k-d','LineWidth',1.5,'MarkerSize',8)

axis([10 10^7 10^(-0) 10^2])

set(gca,'FontSize',fsize)
set(gcf, 'Color', 'w');

 [hleg1, hobj1]=legend('$C_f=10^{-6}$' ,'$C_f=10^{-4}$',...
     '$C_f=10^{-2}$','$C_f=1$');
 h=legend;
 set(h, 'interpreter', 'latex','Fontsize',fsize);
 set(hleg1,'units','normalized','position',[0 0 .3 .1])
 set(gca,'XTick',[10 10^3 10^5 10^7])
 
 xlabel('DOF','interpreter', 'latex','Fontsize',fsize);
 ylabel('estimator/error','interpreter', 'latex','Fontsize',fsize);
 
 
 
 
 
%  load testcf_1_eps1.mat
% N_1=errinfo(:,1);
% err_1=errinfo(:,2);
% est_1=errinfo(:,3)+errinfo(:,4)+errinfo(:,5);
% eff_inf_1=max(err_1./est_1)
% 
% 
% % load testcf_.01.mat
% % N_01=errinfo(:,1);
% % lg=length(N_01);
% % N_01=N_01(1:lg-1);
% % err_01=errinfo(1:lg-1,2);
% % est_01=errinfo(1:lg-1,3)+errinfo(1:lg-1,4)+errinfo(1:lg-1,5);
% % eff_inf_01=max(err_01./est_01)
% % 
% % 
% % 
% % 
% % load testcf_.0001.mat
% % N_0001=errinfo(:,1);
% % err_0001=errinfo(:,2);
% % est_0001=errinfo(:,3)+errinfo(:,4)+errinfo(:,5);
% % eff_inf_0001=max(err_0001./est_0001)
% 
% 
% 
% load testcf_1e-6_eps1.mat
% N_000001=errinfo(:,1);
% err_000001=errinfo(:,2);
% est_000001=errinfo(:,3)+errinfo(:,4)+errinfo(:,5);
% 
% eff_inf_000001=max(err_000001./est_000001)
% 
% % load testcf_.00000001.mat
% % N_00000001=errinfo(:,1);
% % err_00000001=errinfo(:,2);
% % est_00000001=errinfo(:,3)+errinfo(:,4)+errinfo(:,5);
% % 
% % eff_inf_00000001=max(err_00000001./est_00000001)
% 
% 
% hFig=figure(3)
% set(hFig, 'units','normalized','Position', pstn)
% 
% loglog(N_000001, err_000001,'k-+','LineWidth',1.5,'MarkerSize',8)
% hold on
% %loglog(N_0001, err_0001,'k-o','LineWidth',1.5,'MarkerSize',8)
% loglog(N_1, 150*log(N_1)./N_1,'k--','LineWidth',1.5, 'MarkerSize',8)
% %loglog(N_01, err_01,'k-*','LineWidth',1.5,'MarkerSize',8)
% loglog(N_1, err_1,'k-d','LineWidth',1.5,'MarkerSize',8)
% %loglog(N_1, est_1, 'k-x','LineWidth',1.5,'MarkerSize',8)
% 
% 
% %loglog(N_00000001,err_00000001,'k-v','LineWidth',1.5,'MarkerSize',8)
% axis([10 10^7 10^(-4) 10])
% 
% set(gca,'FontSize',fsize)
% set(gcf, 'Color', 'w');
% 
%  [hleg1, hobj1]=legend('$C_f=10^{-6}$' ,'$C_f=10^{-4}$',...
%      '$\log(DOF)/DOF$','$C_f=10^{-2}$','$C_f=1$','est, $C_f=1$');
%  h=legend;
%  set(h, 'interpreter', 'latex','Fontsize',fsize);
%  set(hleg1,'units','normalized','position',[0 0 .39 .1])
%  set(gca,'XTick',[10 10^3 10^5 10^7])
%  
%  xlabel('DOF','interpreter', 'latex','Fontsize',fsize);
%  ylabel('error','interpreter', 'latex','Fontsize',fsize);

 
 load fdiscont_1e-0_newer.mat
 
 N_1=errinfo(:,1);
 residest_1=errinfo(:,3);
 qest_1=errinfo(:,4);
 qmest_1=errinfo(:,5);
 sumqest_1=errinfo(:,9);
 sumqmest_1=errinfo(:,10);
 maxqest_1=errinfo(:,11);
 maxqmest_1=errinfo(:,12);
 
  load fdiscont_1e-4_newer.mat

 N_0001=errinfo(:,1);
 residest_0001=errinfo(:,3);
 qest_0001=errinfo(:,4);
 qmest_0001=errinfo(:,5);
 sumqest_0001=errinfo(:,9);
 sumqmest_0001=errinfo(:,10);
 maxqest_0001=errinfo(:,11);
 maxqmest_0001=errinfo(:,12);
 
 
 hFig=figure(4)
set(hFig, 'units','normalized','Position', pstn)

loglog(N_0001, qest_0001,'k-*','LineWidth',1.5,'MarkerSize',8)
hold on
%loglog(N_0001, err_0001,'k-o','LineWidth',1.5,'MarkerSize',8)
loglog(N_0001, 3000*log(N_0001)./N_0001,'k--','LineWidth',1.5, 'MarkerSize',8)
loglog(N_0001, residest_0001,'k-+','LineWidth',1.5,'MarkerSize',8)
loglog(N_0001, qmest_0001,'k-d','LineWidth',1.5,'MarkerSize',8)
%loglog(N_1, est_1, 'k-x','LineWidth',1.5,'MarkerSize',8)


%loglog(N_00000001,err_00000001,'k-v','LineWidth',1.5,'MarkerSize',8)
axis([10 10^7 10^(-4) 10])

set(gca,'FontSize',fsize)
set(gcf, 'Color', 'w');

 [hleg1, hobj1]=legend('$\eta^q$',...
     '$\log(DOF)/DOF$','$\eta^\infty$','$\eta^{q-1}$');
 h=legend;
 set(h, 'interpreter', 'latex','Fontsize',fsize);
 set(hleg1,'units','normalized','position',[0 0 .42 .1])
 set(gca,'XTick',[10 10^3 10^5 10^7])
 
 xlabel('DOF','interpreter', 'latex','Fontsize',fsize);
 ylabel('estimator','interpreter', 'latex','Fontsize',fsize);
 
 
  hFig=figure(5)
set(hFig, 'units','normalized','Position', pstn)

loglog(N_0001, sumqest_0001,'k-+','LineWidth',1.5,'MarkerSize',8)
hold on
%loglog(N_0001, err_0001,'k-o','LineWidth',1.5,'MarkerSize',8)
loglog(N_0001, 3000*log(N_0001)./N_0001,'k--','LineWidth',1.5, 'MarkerSize',8)
loglog(N_0001, qest_0001,'k-*','LineWidth',1.5,'MarkerSize',8)
loglog(N_0001, maxqest_0001,'k-d','LineWidth',1.5,'MarkerSize',8)
%loglog(N_1, est_1, 'k-x','LineWidth',1.5,'MarkerSize',8)


%loglog(N_00000001,err_00000001,'k-v','LineWidth',1.5,'MarkerSize',8)
axis([10 10^7 10^(-2) 10^2])

set(gca,'FontSize',fsize)
set(gcf, 'Color', 'w');

%For publication
[hleg1, hobj1]=legend('$\eta^q_{\Sigma}$',...
     '$\log(DOF)/DOF$','$\eta^q$','$\|\mu^q\|_{\infty ;\Omega}$');

h=legend;
 set(h, 'interpreter', 'latex','Fontsize',fsize);
 set(hleg1,'units','normalized','position',[0 0 .42 .1])
 set(gca,'XTick',[10 10^3 10^5 10^7])
 
 xlabel('DOF','interpreter', 'latex','Fontsize',fsize);
 ylabel('estimator','interpreter', 'latex','Fontsize',fsize);
 
 hFig=figure(6)
set(hFig, 'units','normalized','Position', pstn)

loglog(N_1, qest_1,'k-*','LineWidth',1.5,'MarkerSize',8)
hold on
%loglog(N_0001, err_0001,'k-o','LineWidth',1.5,'MarkerSize',8)
loglog(N_1, 1.50*log(N_1)./N_1,'k--','LineWidth',1.5, 'MarkerSize',8)
loglog(N_1, residest_1,'k-+','LineWidth',1.5,'MarkerSize',8)
loglog(N_1, qmest_1,'k-d','LineWidth',1.5,'MarkerSize',8)
%loglog(N_1, est_1, 'k-x','LineWidth',1.5,'MarkerSize',8)


%loglog(N_00000001,err_00000001,'k-v','LineWidth',1.5,'MarkerSize',8)
axis([10 10^7 10^(-7) .1])

set(gca,'FontSize',fsize)
set(gcf, 'Color', 'w');

 [hleg1, hobj1]=legend('$\eta^q$',...
     '$\log(DOF)/DOF$','$\eta^\infty$','$\eta^{q-1}$');
 h=legend;
 set(h, 'interpreter', 'latex','Fontsize',fsize);
 set(hleg1,'units','normalized','position',[0 0 .39 .1])
 set(gca,'XTick',[10 10^3 10^5 10^7])
 
 xlabel('DOF','interpreter', 'latex','Fontsize',fsize);
 ylabel('estimator','interpreter', 'latex','Fontsize',fsize);
 
  hFig=figure(7)
set(hFig, 'units','normalized','Position', pstn)

loglog(N_1, maxqest_1,'k-d','LineWidth',1.5,'MarkerSize',8)
hold on
%loglog(N_0001, err_0001,'k-o','LineWidth',1.5,'MarkerSize',8)
loglog(N_1, 1.5*log(N_1)./N_1,'k--','LineWidth',1.5, 'MarkerSize',8)
loglog(N_1, qest_1,'k-*','LineWidth',1.5,'MarkerSize',8)
loglog(N_1, sumqest_1,'k-+','LineWidth',1.5,'MarkerSize',8)
%loglog(N_1, est_1, 'k-x','LineWidth',1.5,'MarkerSize',8)


%loglog(N_00000001,err_00000001,'k-v','LineWidth',1.5,'MarkerSize',8)
axis([10 10^7 10^(-5) 1])

set(gca,'FontSize',fsize)
set(gcf, 'Color', 'w');

 [hleg1, hobj1]=legend('$\|\mu^q\|_{\infty ;\Omega}$', ...
     '$\log(DOF)/DOF$','$\eta^q$','$\eta^q_{\Sigma}$');
 h=legend;
 set(h, 'interpreter', 'latex','Fontsize',fsize);
 set(hleg1,'units','normalized','position',[0 0 .42 .1])
 set(gca,'XTick',[10 10^3 10^5 10^7])
 
 xlabel('DOF','interpreter', 'latex','Fontsize',fsize);
 ylabel('estimator','interpreter', 'latex','Fontsize',fsize);
 
 load pb_eps1e-6_new.mat
 
 N_pb=errinfo(:,1);
 err_pb=errinfo(:,2);
 res_pb=errinfo(:,3);
 quadest_pb=errinfo(:,4)+errinfo(:,5);
 
  hFig=figure(8)
set(hFig, 'units','normalized','Position', pstn)

loglog(N_pb, err_pb,'k-d','LineWidth',1.5,'MarkerSize',8)
hold on
%loglog(N_0001, err_0001,'k-o','LineWidth',1.5,'MarkerSize',8)
loglog(N_pb, 15*log(N_pb)./N_pb,'k--','LineWidth',1.5, 'MarkerSize',8)
loglog(N_pb, res_pb,'k-*','LineWidth',1.5,'MarkerSize',8)
loglog(N_pb, quadest_pb,'k-+','LineWidth',1.5,'MarkerSize',8)
%loglog(N_1, est_1, 'k-x','LineWidth',1.5,'MarkerSize',8)


%loglog(N_00000001,err_00000001,'k-v','LineWidth',1.5,'MarkerSize',8)
axis([10 10^7 10^(-5) 10])

set(gca,'FontSize',fsize)
set(gcf, 'Color', 'w');

 [hleg1, hobj1]=legend('$\|u-u_h\|_{\infty ;\Omega}$', ...
     '$\log(DOF)/DOF$','$\eta^\infty$','$\eta^q+\eta^{q-1}$');
 h=legend;
 set(h, 'interpreter', 'latex','Fontsize',fsize);
 set(hleg1,'units','normalized','position',[0 0 .42 .1])
 set(gca,'XTick',[10 10^3 10^5 10^7])
 
 xlabel('DOF','interpreter', 'latex','Fontsize',fsize);
 ylabel('estimator','interpreter', 'latex','Fontsize',fsize);
 
 
  load pb_lshape_2_new.mat
N_pb=errinfo(:,1);
 err_pb=errinfo(:,2);
 res_pb=errinfo(:,3);
 quadest_pb=errinfo(:,4)+errinfo(:,5);
 
  hFig=figure(9)
set(hFig, 'units','normalized','Position', pstn)

%loglog(N_pb, err_pb,'k-d','LineWidth',1.5,'MarkerSize',8)
loglog(N_pb, res_pb,'k-*','LineWidth',1.5,'MarkerSize',8)
hold on
%loglog(N_0001, err_0001,'k-o','LineWidth',1.5,'MarkerSize',8)
loglog(N_pb, 150*log(N_pb)./N_pb,'k--','LineWidth',1.5, 'MarkerSize',8)
loglog(N_pb, quadest_pb,'k-+','LineWidth',1.5,'MarkerSize',8)
%loglog(N_1, est_1, 'k-x','LineWidth',1.5,'MarkerSize',8)


%loglog(N_00000001,err_00000001,'k-v','LineWidth',1.5,'MarkerSize',8)
axis([10 10^7 10^(-5) 10])

set(gca,'FontSize',fsize)
set(gcf, 'Color', 'w');

 [hleg1, hobj1]=legend('$\eta^\infty$', ...
     '$\log(DOF)/DOF$','$\eta^q+\eta^{q-1}$');
 h=legend;
 set(h, 'interpreter', 'latex','Fontsize',fsize);
 set(hleg1,'units','normalized','position',[0 0 .42 .1])
 set(gca,'XTick',[10 10^3 10^5 10^7])
 
 xlabel('DOF','interpreter', 'latex','Fontsize',fsize);
 ylabel('estimator','interpreter', 'latex','Fontsize',fsize);
 
 
% load testcf_nonlin_1.mat
% N_1=errinfo(:,1);
% err_1=errinfo(:,2);
% est_1=errinfo(:,3)+errinfo(:,4)+errinfo(:,5);
% eff_inf_1=max(err_1./est_1)
% 
% 
% load testcf_nonlin_1e-2.mat
% N_01=errinfo(:,1);
% lg=length(N_01);
% N_01=N_01(1:lg-1);
% err_01=errinfo(1:lg-1,2);
% est_01=errinfo(1:lg-1,3)+errinfo(1:lg-1,4)+errinfo(1:lg-1,5);
% eff_inf_01=max(err_01./est_01)
% 
% load testcf_nonlin_1e-4.mat
% N_0001=errinfo(:,1);
% err_0001=errinfo(:,2);
% est_0001=errinfo(:,3)+errinfo(:,4)+errinfo(:,5);
% eff_inf_0001=max(err_0001./est_0001)
% 
% 
% 
% load testcf_.000001.mat
% N_000001=errinfo(:,1);
% err_000001=errinfo(:,2);
% est_000001=errinfo(:,3)+errinfo(:,4)+errinfo(:,5);
% 
% eff_inf_000001=max(err_000001./est_000001)
% 
% load testcf_.00000001.mat
% N_00000001=errinfo(:,1);
% err_00000001=errinfo(:,2);
% est_00000001=errinfo(:,3)+errinfo(:,4)+errinfo(:,5);
% 
% eff_inf_00000001=max(err_00000001./est_00000001)
% 
% 
% hFig=figure(3)
% set(hFig, 'units','normalized','Position', pstn)
% 
% loglog(N_1, err_1,'k-d','LineWidth',1.5,'MarkerSize',8)
% hold on
% loglog(N_01, err_01,'k-*','LineWidth',1.5,'MarkerSize',8)
% loglog(N_1, 150*log(N_1)./N_1,'k--','LineWidth',1.5, 'MarkerSize',8)
% loglog(N_0001, err_0001,'k-o','LineWidth',1.5,'MarkerSize',8)
% loglog(N_000001, err_000001,'k-+','LineWidth',1.5,'MarkerSize',8)
% loglog(N_00000001,err_00000001,'k-v','LineWidth',1.5,'MarkerSize',8)
% axis([10 10^7 10^(-4) 10])
% 
% set(gca,'FontSize',fsize)
% set(gcf, 'Color', 'w');
% 
%  [hleg1, hobj1]=legend('$C_f=1$' ,'$C_f=10^{-2}$',...
%      '$\log(DOF)/DOF$','$C_f=10^{-4}$','$C_f=10^{-6}$');
%  h=legend;
%  set(h, 'interpreter', 'latex','Fontsize',fsize);
%  set(hleg1,'units','normalized','position',[0 0 .42 .1])
%  
%  
%  xlabel('DOF','interpreter', 'latex','Fontsize',fsize);
%  ylabel('error','interpreter', 'latex','Fontsize',fsize);



 
 
% hFig = figure(3);
% set(hFig, 'units','normalized','Position', pstn)
% 
% 
% loglog(h, 2*h.^2, 'k--','LineWidth',1.5,'MarkerSize',6)
% hold on
% loglog(h, q0, 'k-o','LineWidth',1.5,'MarkerSize',6)
% loglog(h, q1, 'k-x','LineWidth',1.5,'MarkerSize',6)
% loglog(h, .5*h.^3, 'k:', 'LineWidth',1.5,'MarkerSize',6)
% loglog(h, q2, 'k-*', 'LineWidth',1.5,'MarkerSize',6)
% set(gca,'FontSize',fsize)
% set(gcf, 'Color', 'w');
% 
% [hleg1, hobj1]=legend('$h^2$','$\|q-q_h\|_{L_2},$ $f_h=f^\ell$',  '$\|q-q_h\|_{L_2},$ $f_h=\mu_h f^\ell$',...
%  '$h^3$','$\|q-q_h\|_{L_2}$ parametric' );
% h=legend;
% set(h, 'interpreter', 'latex','Fontsize',fsize);
% set(hleg1,'units','normalized','position',[0 0 .42 .1])
% 
% 
% xlabel('$h$','interpreter', 'latex','Fontsize',fsize);
% ylabel('error','interpreter', 'latex','Fontsize',fsize);