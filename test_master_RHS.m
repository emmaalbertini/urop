

%% Set-up

close all; clear; clc;

set(0,'DefaultLineLineWidth',2);
set(0,'defaultAxesFontSize',24);

E_yellow = [224 199 5]/255;
M_orange = [224 135 0]/255;
M_red = [255,0,0]/255;
lightgrey = [217 217 217]/255;
mediumgrey = [140 140 140]/255;
darkgrey = [25 25 25]/255;
R_blue = [0 114 189]/255;
G_green = [119 172 48]/255;
H_yellow = [237 177 32]/255;

tic % Start timer

%%%%% IDEA: Can I save variables and ODE values as pairs? I.e. one ODE
%%%%% equation for R, however each variable has 2 slots, one for E and one
%%%%% for M. And then extend this to N*cell types...
%% Import

clc; close all; clear all;

% Master RHS
Ecell_config;
Mrcell_config;
Mpcell_config;
Mrpcell_config;

%% Parameters

% Parameters
p.mass = 10^8;
p.N = 10^4;
p.zr = 0.01;
p.zp = 0.01;
p.s = 10;
p.v_e = 1000;
p.K_e = 100;
p.n_R = 22377;
p.n_C = 900;
p.n_P = 24804;
p.n_Q = 900;
p.kb_TX = 10;
p.ku_TX = 1;
p.v_TX = 216000;
p.K_TX = 500;
p.K_Q = 150000;
p.hQ = 4;
p.m_deg = 12;
p.kb_TL = 10;
p.ku_TL = 1;
p.v_TL = 72000;
p.K_TL = 400;
p.plasmid = 1;
p.n_H = 900;
p.prom_plus = 100;
p.prom_minus = 1;
p.RBS_plus = 100;
p.RBS_minus = 1;
p.Mfactor = 0;

bioreactor.Ecell = Ecell;
bioreactor.Mrcell = Mrcell;
bioreactor.Mpcell = Mpcell;
bioreactor.Mcell = Mcell;

bioreactor.y_order = {'Ecell','Mrcell','Mpcell','Mcell'};


%% Initial conditions

Ecell_0 = p.N;
Mrcell_0 = p.N - Ecell_0;
Mpcell_0 = p.N - Ecell_0;
Mcell_0 = p.N - Ecell_0;
e_0 = 0; % Energy
TX_0 = [1,1,1,1,1]; % RNApol-DNA complex
m_0 = [1,1,1,1,1]; % mRNAs
TL_0 = [1,1,1,1,1]; % Ribosome-mRNA complex
p_0 = [30, 6400, 1, 8800, 6400]; % Proteins

var = [e_0, TX_0(1:4), m_0(1:4), TL_0(1:4), p_0(1:4),...
         TX_0(end), m_0(end), TL_0(end), p_0(end)];     
var_E = [Ecell_0, var];
var_Mr = [Mrcell_0, var];
var_Mp = [Mpcell_0, var];
var_M = [Mcell_0, var];

x0 = [var_E, var_Mr, var_Mp, var_M];

%% ODE

% tspan = 1:10:50;
tspan = [0,250];

[T,Y] = ode15s(@(t,y) master_RHS(t, y, p, bioreactor),...
               tspan,...
               x0);

%% EM plot

close all
set(gcf, 'Position',  [0, 100, 800, 640])

plot(T, Y(:,1))
hold on
plot(T, Y(:,23))
plot(T, Y(:,45))
plot(T, Y(:,67))
legend('E','Mr','Mp', 'M')
xlabel('Time / h')
ylabel('Number of cells')
grid on
axis square
%% H plot

close all
set(gcf, 'Position',  [0, 100, 800, 640])

plot(T, Y(:,22))
hold on
plot(T, Y(:,44))
plot(T, Y(:,66))
plot(T, Y(:,88))
legend('H_E','H_{Mr}','H_{Mp}', 'H_M')
xlabel('Time / h')
ylabel('Number of cells')
grid on
axis square
%% m_H plot

close all
set(gcf, 'Position',  [0, 100, 800, 640])

plot(T, Y(:,20))
hold on
plot(T, Y(:,42))
plot(T, Y(:,64))
plot(T, Y(:,86))
xlim([0,50])
legend('mH_E','mH_{Mr}','mH_{Mp}', 'mH_M')
xlabel('Time / h')
ylabel('Number of cells')
grid on
axis square


