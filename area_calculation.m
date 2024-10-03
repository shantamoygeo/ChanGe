clc
clear all
data = xlsread('Section.xlsx'); %add your section data here
[A, Q, D, V, W, h] = crosssection_area(data, 'slope',0.008, 'manning', 0.03);
