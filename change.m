function [Q,V,d,WP,Rh] = change(data,m,S)

% This function calculates several hydraulic and hydrological parameters 
% for natura channel. This is an updated version with signficant 
% simplification of the Change.
%
% Guha, S., Singh, A. and Kaushal, R.K., 2020. ChanGe: a MATLAB-based tool 
% for calculation of channel hydrological parameters. Current Science
% (00113891), 119(5).
%
%   Inputs
%
%       Data - Two columns input needed. First column is distance in X axis
%       and the second column is the evelation.
%
%       S - Slope, m - Manning's roughness 
%
%   Output
%       
%       Q - Discharge
%
%       V - Velocity
%       d - Average depth
%       WP - Wetted perimeter
%       Rh - Hydraulic radius
%
%
% Example
%   
%       data = xlsread('data.xlsx');
%       m = 0.045;
%       S = 0.003;
%       [Q,V,d,WP,Rh] = change(data,m,S);



%% Part 1: Assigning X and Y values for the data and creating the cross section

X=data(:,1);
Y=data(:,2);

figure;

plot(X,Y,'o-','LineWidth',0.5)
grid on
xlabel('X');
ylabel('Y');
title('River Cross-Section');


%% Part 2: Seclection of a point by user and finding the neares point

disp('Click on a point near the any point at one bank');
[x_user,y_user] = ginput(1);  % User clicks on figure

dists=sqrt((X-x_user).^2+(Y-y_user).^2); % Distance from each point on cross section
[~,nearest_idx]=min(dists);  % Index of nearest point of the click

hold on;
plot(X(nearest_idx), Y(nearest_idx),'ro','MarkerSize',10,'LineWidth',2); % Marking the nearest point


%% Part 4: Horizonal water level
y_line=Y(nearest_idx); % Y of selected point
x_min=min(X);
x_max=max(X);

hold on;
% Plotting selected point
plot([x_min x_max],[y_line y_line],'b--','LineWidth',1.5);


%% Part 5: Finding out the intersection point between the Yline and cross section
X_dom=[min(X) max(X)];  % span entire X
Y_dom=[y_line y_line];  % Y is fixed

% Finding the intersection points
[xi,yi]=polyxpoly(X, Y, X_dom, Y_dom);
idx_within=(X>=xi(1)) & (X<=xi(2)) & (Y<=yi(1));
X_cross=X(idx_within);
Y_cross=Y(idx_within);
inter_point=horzcat(xi, yi);


% Finding the enique coordinate and concatenating to the cross section
cross=[X_cross Y_cross];
row_unique=setdiff(inter_point, cross,'rows');

if cross(1,2)==row_unique(1,2)
    cross_new=vertcat(cross,row_unique);
else 
    cross_new=vertcat(row_unique,cross);
end

% Extracting columns for further calculations
X_cross_fin=cross_new(:,1);
Y_cross_fin=cross_new(:,2);

% Plot the cross-section
plot(X_cross_fin,Y_cross_fin,'b-o','LineWidth',1.5,'MarkerSize',5);
title('Extracted Cross-Section');

%% Calculation for all the hydraulic parameters

% Calculate cross-section area
area_box=(max(X_cross_fin)-min(X_cross_fin))*max(Y_cross_fin);
area_under_curve=trapz(X_cross_fin,Y_cross_fin);
A = area_box-area_under_curve;

% Wetted perimeter
dx=diff(X_cross_fin);     % difference in X between consecutive points
dy=diff(Y_cross_fin);     % difference in Y between consecutive points
segment_distances=sqrt(dx.^2+dy.^2);
WP = sum(segment_distances);

% Average depth
d=mean(y_line-Y_cross_fin);

% Hydraulic radius
Rh=A/WP;

% Discharge calculation
V=(1/m)*((Rh)^(2/3))*(S)^0.5;
Q=A*V;
