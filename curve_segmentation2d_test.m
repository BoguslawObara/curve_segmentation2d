%% clear
clc; clear all; close all;

%% path
addpath('./lib/');

%% generate curve
m = curveCircle2('close');

%% segment curve 
s = curve_segmentation2d(m);

%% Plot
figure;
plot(m(:,1), m(:,2),'*'); hold on;
plot(m(s(:),1), m(s(:),2),'ro','LineWidth',2,'MarkerEdgeColor','r','MarkerSize',10); 