clear; clc;
load PID_Quad_simulations.mat
format long

t           = out.Desired.time;
x_d         = out.Desired.signals.values(:,1);
y_d         = out.Desired.signals.values(:,2);
z_d         = out.Desired.signals.values(:,3);
psi_d       = rad2deg(out.Desired.signals.values(:,4));

%%% Model Based Controller Signals
phi_d_MBC   = rad2deg(out.MBC.signals.values(:,1));
theta_d_MBC = rad2deg(out.MBC.signals.values(:,2));
x_MBC       = out.MBC.signals.values(:,4);
y_MBC       = out.MBC.signals.values(:,5);
z_MBC       = out.MBC.signals.values(:,6);
phi_MBC     = rad2deg(out.MBC.signals.values(:,7));
theta_MBC   = rad2deg(out.MBC.signals.values(:,8));
psi_MBC     = rad2deg(out.MBC.signals.values(:,9));
F_MBC       = out.MBC.signals.values(:,10);
Tphi_MBC    = out.MBC.signals.values(:,11);
Ttheta_MBC  = out.MBC.signals.values(:,12);
Tpsi_MBC    = out.MBC.signals.values(:,13);

ex_MBC     = x_d - x_MBC;
ey_MBC     = y_d - y_MBC;
ez_MBC     = z_d - z_MBC;
ephi_MBC   = phi_d_MBC - phi_MBC;
etheta_MBC = theta_d_MBC - theta_MBC;
epsi_MBC   = psi_d - psi_MBC;

%%% Adaptive Neural Network Controller Signals 
phi_d_ANNC   = rad2deg(out.ANNC.signals.values(:,1));
theta_d_ANNC = rad2deg(out.ANNC.signals.values(:,2));
x_ANNC       = out.ANNC.signals.values(:,4);
y_ANNC       = out.ANNC.signals.values(:,5);
z_ANNC       = out.ANNC.signals.values(:,6);
phi_ANNC     = rad2deg(out.ANNC.signals.values(:,7));
theta_ANNC   = rad2deg(out.ANNC.signals.values(:,8));
psi_ANNC     = rad2deg(out.ANNC.signals.values(:,9));
F_ANNC       = out.ANNC.signals.values(:,10);
Tphi_ANNC    = out.ANNC.signals.values(:,11);
Ttheta_ANNC  = out.ANNC.signals.values(:,12);
Tpsi_ANNC    = out.ANNC.signals.values(:,13);

ex_ANNC     = x_d - x_ANNC;
ey_ANNC     = y_d - y_ANNC;
ez_ANNC     = z_d - z_ANNC;
ephi_ANNC   = phi_d_ANNC - phi_ANNC;
etheta_ANNC = theta_d_ANNC - theta_ANNC;
epsi_ANNC   = psi_d - psi_ANNC;

%%% Model-free PID Controller Signals Comparison
phi_d_MFPID   = rad2deg(out.MFPID.signals.values(:,1));
theta_d_MFPID = rad2deg(out.MFPID.signals.values(:,2));
x_MFPID       = out.MFPID.signals.values(:,4);
y_MFPID       = out.MFPID.signals.values(:,5);
z_MFPID       = out.MFPID.signals.values(:,6);
phi_MFPID     = rad2deg(out.MFPID.signals.values(:,7));
theta_MFPID   = rad2deg(out.MFPID.signals.values(:,8));
psi_MFPID     = rad2deg(out.MFPID.signals.values(:,9));
F_MFPID       = out.MFPID.signals.values(:,10);
Tphi_MFPID    = out.MFPID.signals.values(:,11);
Ttheta_MFPID  = out.MFPID.signals.values(:,12);
Tpsi_MFPID    = out.MFPID.signals.values(:,13);

ex_MFPID     = x_d - x_MFPID;
ey_MFPID     = y_d - y_MFPID;
ez_MFPID     = z_d - z_MFPID;
ephi_MFPID   = phi_d_MFPID - phi_MFPID;
etheta_MFPID = theta_d_MFPID - theta_MFPID;
epsi_MFPID   = psi_d - psi_MFPID;

%%% Nonlinear PID-Type Controller Signals Comparison
phi_d_NLPIDT   = rad2deg(out.NLPIDT.signals.values(:,1));
theta_d_NLPIDT = rad2deg(out.NLPIDT.signals.values(:,2));
x_NLPIDT       = out.NLPIDT.signals.values(:,4);
y_NLPIDT       = out.NLPIDT.signals.values(:,5);
z_NLPIDT       = out.NLPIDT.signals.values(:,6);
phi_NLPIDT     = rad2deg(out.NLPIDT.signals.values(:,7));
theta_NLPIDT   = rad2deg(out.NLPIDT.signals.values(:,8));
psi_NLPIDT     = rad2deg(out.NLPIDT.signals.values(:,9));
F_NLPIDT       = out.NLPIDT.signals.values(:,10);
Tphi_NLPIDT    = out.NLPIDT.signals.values(:,11);
Ttheta_NLPIDT  = out.NLPIDT.signals.values(:,12);
Tpsi_NLPIDT    = out.NLPIDT.signals.values(:,13);

ex_NLPIDT     = x_d - x_NLPIDT;
ey_NLPIDT     = y_d - y_NLPIDT;
ez_NLPIDT     = z_d - z_NLPIDT;
ephi_NLPIDT   = phi_d_NLPIDT - phi_NLPIDT;
etheta_NLPIDT = theta_d_NLPIDT - theta_NLPIDT;
epsi_NLPIDT   = psi_d - psi_NLPIDT;

%%% Nonlinear PID-regulation Controller Signals Comparison
phi_d_NLPIDR   = rad2deg(out.NLPIDR.signals.values(:,1));
theta_d_NLPIDR = rad2deg(out.NLPIDR.signals.values(:,2));
x_NLPIDR       = out.NLPIDR.signals.values(:,4);
y_NLPIDR       = out.NLPIDR.signals.values(:,5);
z_NLPIDR       = out.NLPIDR.signals.values(:,6);
phi_NLPIDR     = rad2deg(out.NLPIDR.signals.values(:,7));
theta_NLPIDR   = rad2deg(out.NLPIDR.signals.values(:,8));
psi_NLPIDR     = rad2deg(out.NLPIDR.signals.values(:,9));
F_NLPIDR       = out.NLPIDR.signals.values(:,10);
Tphi_NLPIDR    = out.NLPIDR.signals.values(:,11);
Ttheta_NLPIDR  = out.NLPIDR.signals.values(:,12);
Tpsi_NLPIDR    = out.NLPIDR.signals.values(:,13);

ex_NLPIDR     = x_d - x_NLPIDR;
ey_NLPIDR     = y_d - y_NLPIDR;
ez_NLPIDR     = z_d - z_NLPIDR;
ephi_NLPIDR   = phi_d_NLPIDR - phi_NLPIDR;
etheta_NLPIDR = theta_d_NLPIDR - theta_NLPIDR;
epsi_NLPIDR   = psi_d - psi_NLPIDR;

%% RMS Computation
%-----------------
interval1 = t>=20 & t<40;   % Undisturbed
interval2 = t>=40 & t<=60; % Dirturbances
%----- x coordinate
rms_ex_MB_1  = rms(ex_MBC(interval1));
rms_ex_ANN_1 = rms(ex_ANNC(interval1));
rms_ex_MF_1  = rms(ex_MFPID(interval1));
rms_ex_NL_1  = rms(ex_NLPIDT(interval1));
rms_ex_NLR_1 = rms(ex_NLPIDR(interval1));
rms_ex_MB_2  = rms(ex_MBC(interval2));
rms_ex_ANN_2 = rms(ex_ANNC(interval2));
rms_ex_MF_2  = rms(ex_MFPID(interval2));
rms_ex_NL_2  = rms(ex_NLPIDT(interval2));
rms_ex_NLR_2 = rms(ex_NLPIDR(interval2));
%----- y coordinate
rms_ey_MB_1  = rms(ey_MBC(interval1));
rms_ey_ANN_1 = rms(ey_ANNC(interval1));
rms_ey_MF_1  = rms(ey_MFPID(interval1));
rms_ey_NL_1  = rms(ey_NLPIDT(interval1));
rms_ey_NLR_1 = rms(ey_NLPIDR(interval1));
rms_ey_MB_2  = rms(ey_MBC(interval2));
rms_ey_ANN_2 = rms(ey_ANNC(interval2));
rms_ey_MF_2  = rms(ey_MFPID(interval2));
rms_ey_NL_2  = rms(ey_NLPIDT(interval2));
rms_ey_NLR_2 = rms(ey_NLPIDR(interval2));
%----- z coordinate
rms_ez_MB_1  = rms(ez_MBC(interval1));
rms_ez_ANN_1 = rms(ez_ANNC(interval1));
rms_ez_MF_1  = rms(ez_MFPID(interval1));
rms_ez_NL_1  = rms(ez_NLPIDT(interval1));
rms_ez_NLR_1 = rms(ez_NLPIDR(interval1));
rms_ez_MB_2  = rms(ez_MBC(interval2));
rms_ez_ANN_2 = rms(ez_ANNC(interval2));
rms_ez_MF_2  = rms(ez_MFPID(interval2));
rms_ez_NL_2  = rms(ez_NLPIDT(interval2));
rms_ez_NLR_2 = rms(ez_NLPIDR(interval2));
%----- phi angle
rms_ephi_MB_1  = rms(ephi_MBC(interval1));
rms_ephi_ANN_1 = rms(ephi_ANNC(interval1));
rms_ephi_MF_1  = rms(ephi_MFPID(interval1));
rms_ephi_NL_1  = rms(ephi_NLPIDT(interval1));
rms_ephi_NLR_1 = rms(ephi_NLPIDR(interval1));
rms_ephi_MB_2  = rms(ephi_MBC(interval2));
rms_ephi_ANN_2 = rms(ephi_ANNC(interval2));
rms_ephi_MF_2  = rms(ephi_MFPID(interval2));
rms_ephi_NL_2  = rms(ephi_NLPIDT(interval2));
rms_ephi_NLR_2 = rms(ephi_NLPIDR(interval2));
%----- theta angle
rms_etheta_MB_1  = rms(etheta_MBC(interval1));
rms_etheta_ANN_1 = rms(etheta_ANNC(interval1));
rms_etheta_MF_1  = rms(etheta_MFPID(interval1));
rms_etheta_NL_1  = rms(etheta_NLPIDT(interval1));
rms_etheta_NLR_1 = rms(etheta_NLPIDR(interval1));
rms_etheta_MB_2  = rms(etheta_MBC(interval2));
rms_etheta_ANN_2 = rms(etheta_ANNC(interval2));
rms_etheta_MF_2  = rms(etheta_MFPID(interval2));
rms_etheta_NL_2  = rms(etheta_NLPIDT(interval2));
rms_etheta_NLR_2 = rms(etheta_NLPIDR(interval2));
%---- psi angle
rms_epsi_MB_1  = rms(epsi_MBC(interval1));
rms_epsi_ANN_1 = rms(epsi_ANNC(interval1));
rms_epsi_MF_1  = rms(epsi_MFPID(interval1));
rms_epsi_NL_1  = rms(epsi_NLPIDT(interval1));
rms_epsi_NLR_1 = rms(epsi_NLPIDR(interval1));
rms_epsi_MB_2  = rms(epsi_MBC(interval2));
rms_epsi_ANN_2 = rms(epsi_ANNC(interval2));
rms_epsi_MF_2  = rms(epsi_MFPID(interval2));
rms_epsi_NL_2  = rms(epsi_NLPIDT(interval2));
rms_epsi_NLR_2 = rms(epsi_NLPIDR(interval2));
%---- Total thrust
rms_F_MB_1  = rms(F_MBC(interval1));
rms_F_ANN_1 = rms(F_ANNC(interval1));
rms_F_MF_1  = rms(F_MFPID(interval1));
rms_F_NL_1  = rms(F_NLPIDT(interval1));
rms_F_NLR_1 = rms(F_NLPIDR(interval1));
rms_F_MB_2  = rms(F_MBC(interval2));
rms_F_ANN_2 = rms(F_ANNC(interval2));
rms_F_MF_2  = rms(F_MFPID(interval2));
rms_F_NL_2  = rms(F_NLPIDT(interval2));
rms_F_NLR_2 = rms(F_NLPIDR(interval2));
%---- Tau phi
rms_Tphi_MB_1  = rms(Tphi_MBC(interval1));
rms_Tphi_ANN_1 = rms(Tphi_ANNC(interval1));
rms_Tphi_MF_1  = rms(Tphi_MFPID(interval1));
rms_Tphi_NL_1  = rms(Tphi_NLPIDT(interval1));
rms_Tphi_NLR_1 = rms(Tphi_NLPIDR(interval1));
rms_Tphi_MB_2  = rms(Tphi_MBC(interval2));
rms_Tphi_ANN_2 = rms(Tphi_ANNC(interval2));
rms_Tphi_MF_2  = rms(Tphi_MFPID(interval2));
rms_Tphi_NL_2  = rms(Tphi_NLPIDT(interval2));
rms_Tphi_NLR_2 = rms(Tphi_NLPIDR(interval2));
%---- Tau theta
rms_Ttheta_MB_1  = rms(Ttheta_MBC(interval1));
rms_Ttheta_ANN_1 = rms(Ttheta_ANNC(interval1));
rms_Ttheta_MF_1  = rms(Ttheta_MFPID(interval1));
rms_Ttheta_NL_1  = rms(Ttheta_NLPIDT(interval1));
rms_Ttheta_NLR_1 = rms(Ttheta_NLPIDR(interval1));
rms_Ttheta_MB_2  = rms(Ttheta_MBC(interval2));
rms_Ttheta_ANN_2 = rms(Ttheta_ANNC(interval2));
rms_Ttheta_MF_2  = rms(Ttheta_MFPID(interval2));
rms_Ttheta_NL_2  = rms(Ttheta_NLPIDT(interval2));
rms_Ttheta_NLR_2 = rms(Ttheta_NLPIDR(interval2));
%---- Tau psi
rms_Tpsi_MB_1  = rms(Tpsi_MBC(interval1));
rms_Tpsi_ANN_1 = rms(Tpsi_ANNC(interval1));
rms_Tpsi_MF_1  = rms(Tpsi_MFPID(interval1));
rms_Tpsi_NL_1  = rms(Tpsi_NLPIDT(interval1));
rms_Tpsi_NLR_1 = rms(Tpsi_NLPIDR(interval1));
rms_Tpsi_MB_2  = rms(Tpsi_MBC(interval2));
rms_Tpsi_ANN_2 = rms(Tpsi_ANNC(interval2));
rms_Tpsi_MF_2  = rms(Tpsi_MFPID(interval2));
rms_Tpsi_NL_2  = rms(Tpsi_NLPIDT(interval2));
rms_Tpsi_NLR_2 = rms(Tpsi_NLPIDR(interval2));
%------ Creation of vectors
format short
Signal = {'e_x'; 'e_y'; 'e_z'; 'e_phi'; 'e_theta'; 'e_psi'; 'F'; 'Tphi'; 'Ttheta'; 'Tpsi'};
%---- Model-based PID
MB_1  = round([rms_ex_MB_1; rms_ey_MB_1; rms_ez_MB_1; rms_ephi_MB_1; rms_etheta_MB_1; rms_epsi_MB_1; rms_F_MB_1; rms_Tphi_MB_1; rms_Ttheta_MB_1; rms_Tpsi_MB_1],4);
MB_2  = round([rms_ex_MB_2; rms_ey_MB_2; rms_ez_MB_2; rms_ephi_MB_2; rms_etheta_MB_2; rms_epsi_MB_2; rms_F_MB_2; rms_Tphi_MB_2; rms_Ttheta_MB_2; rms_Tpsi_MB_2],4);
%---- ANN + PID
ANN_1 = round([rms_ex_ANN_1; rms_ey_ANN_1; rms_ez_ANN_1; rms_ephi_ANN_1; rms_etheta_ANN_1; rms_epsi_ANN_1; rms_F_ANN_1; rms_Tphi_ANN_1; rms_Ttheta_ANN_1; rms_Tpsi_ANN_1],4);
ANN_2 = round([rms_ex_ANN_2; rms_ey_ANN_2; rms_ez_ANN_2; rms_ephi_ANN_2; rms_etheta_ANN_2; rms_epsi_ANN_2; rms_F_ANN_2; rms_Tphi_ANN_2; rms_Ttheta_ANN_2; rms_Tpsi_ANN_2],4);
%---- Model-free PID
MF_1  = round([rms_ex_MF_1; rms_ey_MF_1; rms_ez_MF_1; rms_ephi_MF_1; rms_etheta_MF_1; rms_epsi_MF_1; rms_F_MF_1; rms_Tphi_MF_1; rms_Ttheta_MF_1; rms_Tpsi_MF_1],4);
MF_2  = round([rms_ex_MF_2; rms_ey_MF_2; rms_ez_MF_2; rms_ephi_MF_2; rms_etheta_MF_2; rms_epsi_MF_2; rms_F_MF_2; rms_Tphi_MF_2; rms_Ttheta_MF_2; rms_Tpsi_MF_2],4);
%---- Nonlinear PID
NL_1  = round([rms_ex_NL_1; rms_ey_NL_1; rms_ez_NL_1; rms_ephi_NL_1; rms_etheta_NL_1; rms_epsi_NL_1; rms_F_NL_1; rms_Tphi_NL_1; rms_Ttheta_NL_1; rms_Tpsi_NL_1],4);
NL_2  = round([rms_ex_NL_2; rms_ey_NL_2; rms_ez_NL_2; rms_ephi_NL_2; rms_etheta_NL_2; rms_epsi_NL_2; rms_F_NL_2; rms_Tphi_NL_2; rms_Ttheta_NL_2; rms_Tpsi_NL_2],4);
%---- Nonlinear PID regulation
NLR_1 = round([rms_ex_NLR_1; rms_ey_NLR_1; rms_ez_NLR_1; rms_ephi_NLR_1; rms_etheta_NLR_1; rms_epsi_NLR_1; rms_F_NLR_1; rms_Tphi_NLR_1; rms_Ttheta_NLR_1; rms_Tpsi_NLR_1],4);
NLR_2 = round([rms_ex_NLR_2; rms_ey_NLR_2; rms_ez_NLR_2; rms_ephi_NLR_2; rms_etheta_NLR_2; rms_epsi_NLR_2; rms_F_NLR_2; rms_Tphi_NLR_2; rms_Ttheta_NLR_2; rms_Tpsi_NLR_2],4);
%---- Table creation
format long
table(Signal, MF_1, NL_1, ANN_1, MB_1, NLR_1, MF_2, NL_2, ANN_2, MB_2, NLR_2)
TABLE_1 = [MF_1 NL_1 ANN_1 MB_1 NLR_1 MF_2 NL_2 ANN_2 MB_2 NLR_2];
%% Figures
fig_width  = 308*2; % 3.2 inches
fig_heigth = fig_width;
xtxtsize = 14;  
ytxtsize = 14;
axessize = 10;
lgndsize = 12;
% dark green color = [0 0.5216 0.1373];
close all
%% Errores
green = [0 0.5216 0.1373];
orange = [1 0.4353 0];
purple = [0.3176 0 1];
f1 = figure('Position',[100 50 fig_width fig_heigth]);
line = 1;
subplot1 = subplot(3,1,1);
e1 = plot(t,0*ex_MBC,'k','linewidth',line);
hold on
e2 = plot(t,ex_MBC,'b','linewidth',line);
e3 = plot(t,ex_ANNC,'r','linewidth',line);
e4 = plot(t,ex_MFPID,'color', green,'linewidth',line);
e5 = plot(t,ex_NLPIDT,'color', orange,'linewidth',line);
e6 = plot(t,ex_NLPIDR,'color', purple,'linewidth',line);
xlim([0 60])
ylim([-0.5 0.5])
grid minor
set(subplot1,...
    'FontName','Times New Roman',...
    'FontSize',axessize,...
    'XTickLabel',{});
ylabel('$e_x$ [m]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);
hleg1 = legend([e4 e5 e3 e2 e6],'MF','NL','ANN','MB','NLR','Location','North','Orientation','horizontal');
set(hleg1,...
    'Interpreter','latex',...
    'NumColumns',5,...
    'FontSize',lgndsize);
legend boxoff;

subplot2 = subplot(3,1,2);
plot(t,0*ey_MBC,'k','linewidth',line);
hold on
plot(t,ey_MBC,'b','linewidth',line);
plot(t,ey_ANNC,'r','linewidth',line);
plot(t,ey_MFPID,'color', green,'linewidth',line);
plot(t,ey_NLPIDT,'color', orange,'linewidth',line);
plot(t,ey_NLPIDR,'color', purple,'linewidth',line);
xlim([0 60])
ylim([-0.5 0.5])
grid minor
set(subplot2,...
    'FontName','Times New Roman',...
    'FontSize',axessize,...
    'XTickLabel',{});
ylabel('$e_y$ [m]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);

subplot3 = subplot(3,1,3);
plot(t,0*ez_MBC,'k','linewidth',line);
hold on
plot(t,ez_MBC,'b','linewidth',line);
plot(t,ez_ANNC,'r','linewidth',line);
plot(t,ez_MFPID,'color', green,'linewidth',line);
plot(t,ez_NLPIDT,'color', orange,'linewidth',line);
plot(t,ez_NLPIDR,'color', purple,'linewidth',line);
xlim([0 60])
ylim([-0.05 0.05])
grid minor
set(subplot3,...
    'FontName','Times New Roman',...
    'FontSize',axessize);
xlabel('Time [s]',...
       'Interpreter','latex',...
       'FontSize',xtxtsize);
ylabel('$e_z$ [m]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);
%% attitide errors
f2 = figure('Position',[100 50 fig_width fig_heigth]);
line = 1;
subplot1 = subplot(3,1,1);
e1 = plot(t,0*ephi_MBC,'k','linewidth',line);
hold on
e2 = plot(t,ephi_MBC,'b','linewidth',line);
e3 = plot(t,ephi_ANNC,'r','linewidth',line);
e4 = plot(t,ephi_MFPID,'color', green,'linewidth',line);
e5 = plot(t,ephi_NLPIDT,'color', orange,'linewidth',line);
e6 = plot(t,ephi_NLPIDR,'color', purple,'linewidth',line);
xlim([0 60])
ylim([-6 6])
grid minor
set(subplot1,...
    'FontName','Times New Roman',...
    'FontSize',axessize,...
    'XTickLabel',{});
ylabel('$e_{\phi}$ [$^\circ$]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);
hleg1 = legend([e4 e5 e3 e2 e6],'MF','NL','ANN','MB','NLR','Location','North','Orientation','horizontal');
set(hleg1,...
    'Interpreter','latex',...
    'NumColumns',5,...
    'FontSize',lgndsize);
legend boxoff;

subplot2 = subplot(3,1,2);
plot(t,0*etheta_MBC,'k','linewidth',line);
hold on
plot(t,etheta_MBC,'b','linewidth',line);
plot(t,etheta_ANNC,'r','linewidth',line);
plot(t,etheta_MFPID,'color', green,'linewidth',line);
plot(t,etheta_NLPIDT,'color', orange,'linewidth',line);
plot(t,etheta_NLPIDR,'color', purple,'linewidth',line);
xlim([0 60])
ylim([-6 6])
grid minor
set(subplot2,...
    'FontName','Times New Roman',...
    'FontSize',axessize,...
    'XTickLabel',{});
ylabel('$e_{\theta}$ [$^\circ$]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);

subplot3 = subplot(3,1,3);
plot(t,0*epsi_MBC,'k','linewidth',line);
hold on
plot(t,epsi_MBC,'b','linewidth',line);
plot(t,epsi_ANNC,'r','linewidth',line);
plot(t,epsi_MFPID,'color', green,'linewidth',line);
plot(t,epsi_NLPIDT,'color', orange,'linewidth',line);
plot(t,epsi_NLPIDR,'color', purple,'linewidth',line);
xlim([0 60])
ylim([-3 3])
grid minor
set(subplot3,...
    'FontName','Times New Roman',...
    'FontSize',axessize);
xlabel('Time [s]',...
       'Interpreter','latex',...
       'FontSize',xtxtsize);
ylabel('$e_{\psi}$ [$^\circ$]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);
   
%% Control Actions
tau_max    = 0.5;
Thrust_max = 20;
f3 = figure('Position',[100 50 fig_width fig_heigth]);
line = 1;
subplot1 = subplot(4,1,1);
e2 = plot(t,F_MBC,'b','linewidth',line);
hold on
e3 = plot(t,F_ANNC,'r','linewidth',line);
e4 = plot(t,F_MFPID,'color', green,'linewidth',line);
e5 = plot(t,F_NLPIDT,'color', orange,'linewidth',line);
e6 = plot(t,F_NLPIDR,'color', purple,'linewidth',line);
xlim([0 60])
ylim([17 Thrust_max])
grid minor
set(subplot1,...
    'FontName','Times New Roman',...
    'FontSize',axessize,...
    'XTickLabel',{});
ylabel('$F$ [N]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);
hleg1 = legend([e4 e5 e3 e2 e6],'MF','NL','ANN','MB','NLR','Location','North','Orientation','horizontal');
set(hleg1,...
    'Interpreter','latex',...
    'NumColumns',5,...
    'FontSize',lgndsize);
legend boxoff;

subplot2 = subplot(4,1,2);
plot(t,Tphi_MBC,'b','linewidth',line);
hold on
plot(t,Tphi_ANNC,'r','linewidth',line);
plot(t,Tphi_MFPID,'color', green,'linewidth',line);
plot(t,Tphi_NLPIDT,'color', orange,'linewidth',line);
plot(t,Tphi_NLPIDR,'color', purple,'linewidth',line);
xlim([0 60])
ylim([-tau_max tau_max])
grid minor
set(subplot2,...
    'FontName','Times New Roman',...
    'FontSize',axessize,...
    'XTickLabel',{});
ylabel('$\tau_{1}$ [Nm]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);

subplot3 = subplot(4,1,3);
plot(t,Ttheta_MBC,'b','linewidth',line);
hold on
plot(t,Ttheta_ANNC,'r','linewidth',line);
plot(t,Ttheta_MFPID,'color', green,'linewidth',line);
plot(t,Ttheta_NLPIDT,'color', orange,'linewidth',line);
plot(t,Ttheta_NLPIDR,'color', purple,'linewidth',line);
xlim([0 60])
ylim([-tau_max tau_max])
grid minor
set(subplot3,...
    'FontName','Times New Roman',...
    'FontSize',axessize,...
    'XTickLabel',{});
ylabel('$\tau_{2}$ [Nm]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);

subplot4 = subplot(4,1,4);
plot(t,Tpsi_MBC,'b','linewidth',line);
hold on
plot(t,Tpsi_ANNC,'r','linewidth',line);
plot(t,Tpsi_MFPID,'color', green,'linewidth',line);
plot(t,Tpsi_NLPIDT,'color', orange,'linewidth',line);
plot(t,Tpsi_NLPIDR,'color', purple,'linewidth',line);
xlim([0 60])
ylim([-tau_max tau_max])
grid minor
set(subplot4,...
    'FontName','Times New Roman',...
    'FontSize',axessize);
xlabel('Time [s]',...
       'Interpreter','latex',...
       'FontSize',xtxtsize);
ylabel('$\tau_{3}$ [Nm]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);
%% 3D Trajectory Tracking
f4 = figure('Position',[100 50 fig_width fig_width]);
axes1 = axes('Parent',f4);
hold(axes1,'on');
plot3(x_d,y_d,z_d,'k','linewidth',1)
hold on
plot3(x_MFPID,y_MFPID,z_MFPID,'color', green,'linewidth',1)
plot3(x_NLPIDT,y_NLPIDT,z_NLPIDT,'color', orange,'linewidth',1)
plot3(x_ANNC,y_ANNC,z_ANNC,'r','linewidth',1)
plot3(x_MBC,y_MBC,z_MBC,'b','linewidth',1)
plot3(x_NLPIDR,y_NLPIDR,z_NLPIDR,'color', purple,'linewidth',1)
plot3(x_d,y_d,z_d,'k','linewidth',1.5)
xlim([-1.2 1.2])
ylim([-1.2 1.2])
zlim([0 1.2])
grid minor

set(axes1,'FontName',...
    'Times New Roman','FontSize',axessize);
xlabel('$x$ [m]',...
    'Interpreter','latex',...
    'FontSize',xtxtsize);
ylabel('$y$ [m]',...
    'Interpreter','latex',...
    'Fontsize',ytxtsize);
zlabel('$z$ [m]',...
    'Interpreter','latex',...
    'Fontsize',ytxtsize);
hleg1 = legend('Desired','MF','NL','ANN','MB','NLR',...
    'Location','Best',...
    'Orientation','horizontal');
set(hleg1,...
    'Interpreter','latex',...
    'NumColumns',1,...
    'FontSize',lgndsize);
legend boxoff
view(axes1,[-45 45]);
%% Position
f5 = figure('Position',[100 50 fig_width fig_heigth]);
line = 1;
subplot1 = subplot(4,1,1);
e1 = plot(t,x_d,'k','linewidth',line);
hold on
e2 = plot(t,x_MBC,'b','linewidth',line);
e3 = plot(t,x_ANNC,'r','linewidth',line);
e4 = plot(t,x_MFPID,'color', green,'linewidth',line);
e5 = plot(t,x_NLPIDT,'color', orange,'linewidth',line);
e6 = plot(t,x_NLPIDR,'color', purple,'linewidth',line);
xlim([0 60])
ylim([-1.2 1.2])
grid minor
set(subplot1,...
    'FontName','Times New Roman',...
    'FontSize',axessize,...
    'XTickLabel',{});
ylabel('$x$ [m]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);
hleg1 = legend([e1 e4 e5 e3 e2 e6],'Desired','MF','NL','ANN','MB','NLR','Location','North','Orientation','horizontal');
set(hleg1,...
    'Interpreter','latex',...
    'NumColumns',6,...
    'FontSize',lgndsize);
legend boxoff;

subplot2 = subplot(4,1,2);
plot(t,y_d,'k','linewidth',line);
hold on
plot(t,y_MBC,'b','linewidth',line);
plot(t,y_ANNC,'r','linewidth',line);
plot(t,y_MFPID,'color', green,'linewidth',line);
plot(t,y_NLPIDT,'color', orange,'linewidth',line);
plot(t,y_NLPIDR,'color', purple,'linewidth',line);
xlim([0 60])
ylim([-1.2 1.2])
grid minor
set(subplot2,...
    'FontName','Times New Roman',...
    'FontSize',axessize,...
    'XTickLabel',{});
ylabel('$y$ [m]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);

subplot3 = subplot(4,1,3);
plot(t,z_d,'k','linewidth',line);
hold on
plot(t,z_MBC,'b','linewidth',line);
plot(t,z_ANNC,'r','linewidth',line);
plot(t,z_MFPID,'color', green,'linewidth',line);
plot(t,z_NLPIDT,'color', orange,'linewidth',line);
plot(t,z_NLPIDR,'color', purple,'linewidth',line);
xlim([0 60])
ylim([0.8 1.05])
grid minor
set(subplot3,...
    'FontName','Times New Roman',...
    'FontSize',axessize,...
    'XTickLabel',{});
ylabel('$z$ [m]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);
   
subplot4 = subplot(4,1,4);
plot(t,psi_d,'k','linewidth',line);
hold on
plot(t,psi_MBC,'b','linewidth',line);
plot(t,psi_ANNC,'r','linewidth',line);
plot(t,psi_MFPID,'color', green,'linewidth',line);
plot(t,psi_NLPIDT,'color', orange,'linewidth',line);
plot(t,psi_NLPIDR,'color', purple,'linewidth',line);
xlim([0 60])
ylim([-3 3])
grid minor
set(subplot4,...
    'FontName','Times New Roman',...
    'FontSize',axessize);
xlabel('Time [s]',...
       'Interpreter','latex',...
       'FontSize',xtxtsize);
ylabel('$\psi$ [$^\circ$]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);
%% Roll and Pitch angles
f6 = figure('Position',[100 50 fig_width fig_heigth]);
line = 1;
%------ Model based
subplot1 = subplot(5,2,7);
desired = plot(t,phi_d_MBC,'k','linewidth',line);
hold on
ann = plot(t,phi_MBC,'r','linewidth',line);
mf = plot(t,phi_MBC,'color', green,'linewidth',line);
nl = plot(t,phi_MBC,'color', orange,'linewidth',line);
nlr = plot(t,phi_MBC,'color', purple,'linewidth',line);
mb = plot(t,phi_MBC,'b','linewidth',line);
xlim([0 60])
lim1 = 15;
ylim([-lim1 lim1])
grid minor
set(subplot1,...
    'FontName','Times New Roman',...
    'FontSize',axessize,...
    'XTickLabel',{});
ylabel('$\phi$ [$^\circ$]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);
hleg1 = legend([desired mf nl ann mb nlr],'Desired','MF','NL','ANN','MB','NLR','Location','North','Orientation','horizontal');
set(hleg1,...
    'Interpreter','latex',...
    'NumColumns',6,...
    'FontSize',lgndsize);
legend boxoff;

subplot2 = subplot(5,2,8);
plot(t,theta_d_MBC,'k','linewidth',line);
hold on
plot(t,theta_MBC,'b','linewidth',line);
xlim([0 60])
ylim([-lim1 lim1])
grid minor
set(subplot2,...
    'FontName','Times New Roman',...
    'FontSize',axessize,...
    'XTickLabel',{});
ylabel('$\theta$ [$^\circ$]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);
%------ ANN
subplot1 = subplot(5,2,5);
plot(t,phi_d_ANNC,'k','linewidth',line);
hold on
plot(t,phi_ANNC,'r','linewidth',line);
xlim([0 60])
ylim([-lim1 lim1])
grid minor
set(subplot1,...
    'FontName','Times New Roman',...
    'FontSize',axessize,...
    'XTickLabel',{});
ylabel('$\phi$ [$^\circ$]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);

subplot2 = subplot(5,2,6);
plot(t,theta_d_ANNC,'k','linewidth',line);
hold on
plot(t,theta_ANNC,'r','linewidth',line);
xlim([0 60])
ylim([-lim1 lim1])
grid minor
set(subplot2,...
    'FontName','Times New Roman',...
    'FontSize',axessize,...
    'XTickLabel',{});
ylabel('$\theta$ [$^\circ$]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);
%------- Model-free PID
subplot1 = subplot(5,2,1);
plot(t,phi_d_MFPID,'k','linewidth',line);
hold on
plot(t,phi_MFPID,'color', green,'linewidth',line);
xlim([0 60])
ylim([-lim1 lim1])
grid minor
set(subplot1,...
    'FontName','Times New Roman',...
    'FontSize',axessize,...
    'XTickLabel',{});
ylabel('$\phi$ [$^\circ$]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);

subplot2 = subplot(5,2,2);
plot(t,theta_d_MFPID,'k','linewidth',line);
hold on
plot(t,theta_MFPID,'color', green,'linewidth',line);
xlim([0 60])
ylim([-lim1 lim1])
grid minor
set(subplot2,...
    'FontName','Times New Roman',...
    'FontSize',axessize,...
    'XTickLabel',{});
ylabel('$\theta$ [$^\circ$]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);
%------ NL PID-Type
subplot1 = subplot(5,2,3);
plot(t,phi_d_NLPIDT,'k','linewidth',line);
hold on
plot(t,phi_NLPIDT,'color', orange,'linewidth',line);
xlim([0 60])
ylim([-lim1 lim1])
grid minor
set(subplot1,...
    'FontName','Times New Roman',...
    'FontSize',axessize,...
    'XTickLabel',{});
ylabel('$\phi$ [$^\circ$]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);

subplot2 = subplot(5,2,4);
plot(t,theta_d_NLPIDT,'k','linewidth',line);
hold on
plot(t,theta_NLPIDT,'color', orange,'linewidth',line);
xlim([0 60])
ylim([-lim1 lim1])
grid minor
set(subplot2,...
    'FontName','Times New Roman',...
    'FontSize',axessize,...
    'XTickLabel',{});
ylabel('$\theta$ [$^\circ$]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);
%------- NL regulation
subplot1 = subplot(5,2,9);
plot(t,phi_d_NLPIDR,'k','linewidth',line);
hold on
plot(t,phi_NLPIDR,'color', purple,'linewidth',line);
xlim([0 60])
ylim([-lim1 lim1])
grid minor
set(subplot1,...
    'FontName','Times New Roman',...
    'FontSize',axessize);
xlabel('Time [s]',...
       'Interpreter','latex',...
       'FontSize',xtxtsize);
ylabel('$\phi$ [$^\circ$]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);

subplot2 = subplot(5,2,10);
plot(t,theta_d_NLPIDR,'k','linewidth',line);
hold on
plot(t,theta_NLPIDR,'color', purple,'linewidth',line);
xlim([0 60])
ylim([-lim1 lim1])
grid minor
set(subplot2,...
    'FontName','Times New Roman',...
    'FontSize',axessize);
xlabel('Time [s]',...
       'Interpreter','latex',...
       'FontSize',xtxtsize);
ylabel('$\theta$ [$^\circ$]',...
       'Interpreter','latex',...
       'FontSize',ytxtsize);
%% Plots
% set(gcf,'renderer','Painters')
% print -f1 -depsc Position_errors_Survey.eps
% print -f2 -depsc Attitude_errors_Survey.eps
% print -f3 -depsc Control_input_Survey.eps
% print -f4 -depsc 3D_Trajectory_Survey.eps
% print -f5 -depsc Position_Survey.eps
% print -f6 -depsc Roll_Pitch_Survey.eps