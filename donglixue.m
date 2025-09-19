clear
clc
tic
%% %%%%%%%%%%%%%%%%%%%%% Large range
load('network.mat'); % load food web
m=length(find(sum(d)==0));
ei = 0.01; % basal species mortality rate
%% establish C-Ctradeoffs among basal species
ci=[0.75;0.25;0.5;1];% basal species colonization rate
H=[0.5 0 0 1;1 0.5 1 1;1 0 0.5 1;0 0 0 0.5]; % competition matrix
%%
cj=0.625; % colonization rate of consumers
ej=0.01; % mortality rate of consumers
u1=0.05; % top-down mortality rate of basal species due to predation
u2=0.05; % top-down mortality rate of consumers when they are fed by other consumers
%% SXR
U1=0:0.01:1; % range of resource productivity
U2=0:0.01:1; % range of ecosystem size
U3=0; % disturbance
HighFF=cell(3,1);
HighFF=equation1(d,U1,U2,U3,ci,cj,ei,ej,H,u1,u2,HighFF);%calling function 'equation1.m'
%% SXD
U1=1; % range of resource productivity
U2=0:0.01:1; % range of ecosystem size
U3=0:0.01:1; % disturbance
HighFF=equation1(d,U1,U2,U3,ci,cj,ei,ej,H,u1,u2,HighFF);%calling function 'equation1.m'
%% DXR
U1=0:0.01:1; % range of resource productivity
U2=1; % range of ecosystem size
U3=0:0.01:1; % disturbance
HighFF=equation1(d,U1,U2,U3,ci,cj,ei,ej,H,u1,u2,HighFF);%calling function 'equation1.m'

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% small range
%% establish C-Ctradeoffs among basal species
ci=linspace(0.45,0.8,m)';% basal species colonization rate
H=[0.5 1 1 1;0 0.5 1 1;0 0 0.5 1;0 0 0 0.5]; % competition matrix
%% SxR
U1=0:0.01:1; % range of resource productivity
U2=0:0.01:1; % range of ecosystem size
U3=0; % disturbance
LowFF=cell(3,1);
LowFF=equation1(d,U1,U2,U3,ci,cj,ei,ej,H,u1,u2,LowFF);%calling function 'equation1.m'
%% SxD
U1=1; % range of resource productivity
U2=0:0.01:1; % range of ecosystem size
U3=0:0.01:1; % disturbance
LowFF=equation1(d,U1,U2,U3,ci,cj,ei,ej,H,u1,u2,LowFF);%calling function 'equation1.m'
%% DxR
U1=0:0.01:1; % range of resource productivity
U2=1; % range of ecosystem size
U3=0:0.01:1; % disturbance
LowFF=equation1(d,U1,U2,U3,ci,cj,ei,ej,H,u1,u2,LowFF);%calling function 'equation1.m'
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting figures
A1=0:0.01:1;%range of x-axis
A2=0:0.01:1;%range of y-axis
sp=5;%maximum food chanin is 5
b=Draw(A1,A2,LowFF,HighFF,shiyitu,wangluo,sp);%calling function "Draw.m" for mapping