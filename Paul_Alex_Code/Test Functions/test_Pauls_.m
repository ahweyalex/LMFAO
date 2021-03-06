% Alexander Moreno
% 06-13-2017
% Testing and creating changes to Paul and Alex's codes
%

clear all; close all; clc;
% example(s)
Ds = 70e-3;     % middle serration depth [m]
M = 50;
er=2.2; %relative permittivity
f0= 433e6; %frequency
%f0 = 18e9;
wm = 5e-3;
we = 20e-3;
sm = 1e-3;
se = 2e-3;
%trial 1
%w = [wm wm wm wm];
%s = [sm sm sm];
% trial 2
%w = [wm wm+we wm+we wm];
%s = [sm sm+we sm];
% trial 3
w = [wm wm wm*1.1];
s = [sm sm];

l= Ds/M;
%h = 0.00254; % 100mil
%h= 0.001524; %substrate height % 60 mil
%h=0.001016; % 40mil
%h=0.000889; % 35mil
%h=0.000762; % 30mil
%h=0.000381; %15mil
%h = 0.000254; %10mil
h=0.000127; % 5mil

Z0 = 50;
ZL = 377;

[Zpaul] = ZParameters(w,h,s,er,f0,l);