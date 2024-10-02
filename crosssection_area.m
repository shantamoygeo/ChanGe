function [cross_section_area,discharge,avg_depth,velocity,wetted_perimeter,hydraulic_radius] = crosssection_area(data,varargin)

% Cross sectional area calculation for the cross section of a natural
% stream using numerical integration method.

% DESCRIPTION
%
% The cross-sectional area calculation for a natural channel is an
% important tool for the calculation of discharge at any point.
% Additionally this tool also calculate the the wetted perimeter and
% velocity when the slope and manning's roughness is given. The main thing
% to remember during the use of this tool is that the points that you click
% on the both bank should be inside channel.
%       Step 1 -  Run the code block given below as per the syntax.
%
%       Step 2 - A figure file will be displayed with the cross profile.
%
%       Step 3 - Click any of the side of the bank. You may click as per
%       the terraces on any of the side. Before you click, please remember
%       that a horizontal line will be generated according to the level of
%       the strath at your chosen side. There should be atleast on point
%       above this straight line on the opposite side of the bank.
%       Otherwise the code will not give proper result.
%
%       Step 4 - Once you click on the both side, the code will calculate
%       series of parameters. These parameters are Cross section area,
%       Wetted perimeter, Hydraulic radius, Flow velocity and Discharge
%
%       Step 5 - Now there will be a figure plotted with the water level
%       and filled with cyan colour. A cursor will be generated. You can
%       click whereever you want and the discharge value will be printed on
%       the graph. You can save the plot.

% Inputs parameter
%
%       Data - This is the main input data for the calculation. This data
%       should contain two columns. The first column should contain all the
%       values along the X axis or distance. The second column should
%       contain the Y axis or the elevation values. Note that the number of
%       X entries and Y entries should be same. These values should be in
%       meters
%
%       % Other Input parameters
%
%       'slope' - This is the slope along the channel. This slope can be
%       calculated from the long profile or from field survey. Note that
%       this slope value should be entered in m/m format
%
%       'manning' - This is the Manning's roughness coefficient. Note that
%       this parameters can be taken from any standard chart.
%
% Output parameters
%       
%       cross_section_area - Calculated cross sectional area.
%
%       discharge - Estimated discharge from Manning's equation. Note that
%       discharge = (1/manning's roughness)*((hydraulic radius)^(2/3))*(slope)^0.5;
%       All the values are in meter.
%
%       avg_depth - Average depth of the channels. The depthe is calculated
%       at each point of the cross section. The depth is calcuated as the
%       distance between the channel and the water level. Therefore the
%       number of the depth values are equal to the number of the point on
%       the cross section-2. This is because the extreme two points are
%       excluded.
%
%
% Example
%   
%       data = xlsread('cross_section.xlsx','sheet1');
%       [A, Q] = crosssection_area(data, 'slope',0.005, 'manning', 0.03);
%
% Process
% The area is calculated using the numerical integration method. The
% integration method gives the area under the curve. Then it is substracted
% from the rectangle covering the area.
%
% Author: Shantamoy Guha (shantamoy.guha@iitgn.ac.in)
% Data: 25. January, 2020


% Parse Inputs
p = inputParser;         
p.FunctionName = 'crosssection_area';
addParameter(p,'slope', @(x) isscalar(x) && x>0);
addParameter(p,'manning', @(x) isscalar(x) && x>0);
parse(p,varargin{:});

%% Calculation of the main cross section
% Inserting the data and dividing into two variable
X = data(:,1);
Y = data(:,2);


% Rounding off
X = round(X,2);
Y = round(Y,2);

% Creating the figure file for pointing the banks

f = figure;
%a = get(f,'Children');
disp('Click on bank')
scatter(X, Y, 'o', 'MarkerFaceColor', 'k')
hold on
plot(X,Y, 'LineWidth', 1.5);
xlabel ('Distance along the bank(m)')
ylabel ('Elevation (m)')
% grid on
a = title('Click near a point on any bank');
[x1,y1] = ginput(1);
line([min(X), max(X)],[y1,y1],'Color','black','LineWidth',1)
a = title('Click on the opposite side where the horizontal black line cuts');
[x2,y2] = ginput(1);
close(f)



% Finding the minimum distance from the selected point left side
if x2>x1
for i = 1:length(X)
    distance_left(i) = sqrt((x1-X(i,1)).^2+(y1-Y(i,1)).^2);
end

minDist_left = min(distance_left);
idx_minDist_left = find(distance_left == minDist_left);


% Finding the minimum distance from the selected point right side
for i = 1:length(X)
    distance_right(i) = sqrt((x2-X(i,1)).^2+(y2-Y(i,1)).^2);
end

minDist_right = min(distance_right);
idx_minDist_right = find(distance_right == minDist_right);

else
    % Finding the minimum distance from the selected point left side
for i = 1:length(X)
    distance_right(i) = sqrt((x1-X(i,1)).^2+(y1-Y(i,1)).^2);
end

minDist_right = min(distance_right);
idx_minDist_right = find(distance_right == minDist_right);


% Finding the minimum distance from the selected point right side
for i = 1:length(X)
    distance_left(i) = sqrt((x2-X(i,1)).^2+(y2-Y(i,1)).^2);
end
minDist_left = min(distance_left);
idx_minDist_left = find(distance_left == minDist_left);  
end



% Adding a data point
if x2>x1 && length(X)~=idx_minDist_right
    X_new = X(idx_minDist_left:idx_minDist_right+1);
    Y_new = Y(idx_minDist_left:idx_minDist_right+1);
    
elseif x2>x1 && length(X)==idx_minDist_right
    X_new = X(idx_minDist_left:idx_minDist_right);
    Y_new = Y(idx_minDist_left:idx_minDist_right);
    
elseif x2<x1 && length(X)~=idx_minDist_right
    X_new = X(idx_minDist_left-1:idx_minDist_right);
    Y_new = Y(idx_minDist_left-1:idx_minDist_right);
    
elseif x2<x1 && length(X)==idx_minDist_right
    X_new = X(idx_minDist_left:idx_minDist_right);
    Y_new = Y(idx_minDist_left:idx_minDist_right);
end


% Finding out the intersection between two curve
X2 = linspace(min(X_new),max(X_new));
if x2>x1
    Y2 = Y_new(1).*ones(size(X2));
else
    Y2 = Y_new(end).*ones(size(X2));
end

C1 = [X_new Y_new]';
C2 = [X2; Y2];

% Finding the intersection points
P = InterX(C1,C2);
P = round(P,2);

if x2>x1
    X_fin = X_new(1:end-1);
    Y_fin = Y_new(1:end-1);
    Area_cal_fin = [X_fin Y_fin];
    extra = P(:,2)';
    Area_cal_fin = vertcat(Area_cal_fin, extra);
else
    X_fin = X_new(2:end);
    Y_fin = Y_new(2:end);
    Area_cal_fin = [X_fin Y_fin];
    extra = P(:,1)';
    Area_cal_fin = vertcat(extra,Area_cal_fin);
end


%% Calculation of area under the curve
Int = trapz(Area_cal_fin(:,1),Area_cal_fin(:,2));
dist = Area_cal_fin(end,1)-Area_cal_fin(1,1);
Area_square = dist*Area_cal_fin (1,2);
cross_section_area = Area_square-Int;


%% Calculation of wetted perimeter, velocity, average depth

% Calculation of the distance between each point.
for i = 1:length(Area_cal_fin)-1
    dist(i) = sqrt((Area_cal_fin(i+1,1)-Area_cal_fin(i,1)).^2+(Area_cal_fin(i+1,2)-Area_cal_fin(i,2)).^2);
end

wetted_perimeter = sum(dist);
hydraulic_radius = round(cross_section_area/wetted_perimeter,2);

slope = p.Results.slope;
manning = p.Results.manning;

velocity = (1/manning)*((hydraulic_radius)^(2/3))*(slope)^0.5;
discharge = cross_section_area*velocity;

dis_for_plot = round (discharge,2);

depths = abs(Area_cal_fin(:,2)-Area_cal_fin(1,2));
avg_depth = sum(depths)/(length(Area_cal_fin)-2);

avg_depth_plot = round (avg_depth,2);

%% Plotting the data

figure
ax1 = axes;
scatter(X, Y, 'o', 'MarkerFaceColor', 'k')
hold on
plot(X,Y, 'LineWidth', 1.5);

hold on

% Colour the paleo cross section
% Make a polygon
ceiling = Area_cal_fin (1,2);  %define top of polygon as top of y axis
xp = [Area_cal_fin(:,1),fliplr(Area_cal_fin(:,1))]; 
yp = [Area_cal_fin(:,2), repmat(ceiling,size(Area_cal_fin(:,1)))]; 
hold(ax1, 'on')
% Fill area above curve
fill(ax1,xp,yp,[0.4 0.8 1])
hold on
xlabel ('Distance along the bank (m)')
ylabel ('Elevation (m)')

[x3,y3] = ginput(1);
str1 = ('Discharge is:');
str2 = num2str(dis_for_plot);
str3 = ('m^3/s');

str = strcat(str1,{'  '},str2,{'  '},str3);

text(x3,y3,str,'FontSize', 15)

[x3,y3] = ginput(1);
str1 = ('Average depth is:');
str2 = num2str(avg_depth_plot);
str3 = ('m');

str = strcat(str1,{'  '},str2,{'  '},str3);

text(x3,y3,str,'FontSize', 15)

end


% A function for calculation of the intersection point
function P = InterX(L1,varargin)
%INTERX Intersection of curves
%   P = INTERX(L1,L2) returns the intersection points of two curves L1 
%   and L2. The curves L1,L2 can be either closed or open and are described
%   by two-row-matrices, where each row contains its x- and y- coordinates.
%   The intersection of groups of curves (e.g. contour lines, multiply 
%   connected regions etc) can also be computed by separating them with a
%   column of NaNs as for example
%
%         L  = [x11 x12 x13 ... NaN x21 x22 x23 ...;
%               y11 y12 y13 ... NaN y21 y22 y23 ...]
%
%   P has the same structure as L1 and L2, and its rows correspond to the
%   x- and y- coordinates of the intersection points of L1 and L2. If no
%   intersections are found, the returned P is empty.
%
%   P = INTERX(L1) returns the self-intersection points of L1. To keep
%   the code simple, the points at which the curve is tangent to itself are
%   not included. P = INTERX(L1,L1) returns all the points of the curve 
%   together with any self-intersection points.
%   
%   Example:
%       t = linspace(0,2*pi);
%       r1 = sin(4*t)+2;  x1 = r1.*cos(t); y1 = r1.*sin(t);
%       r2 = sin(8*t)+2;  x2 = r2.*cos(t); y2 = r2.*sin(t);
%       P = InterX([x1;y1],[x2;y2]);
%       plot(x1,y1,x2,y2,P(1,:),P(2,:),'ro')
%   Author : NS
%   Version: 3.0, 21 Sept. 2010
%   Two words about the algorithm: Most of the code is self-explanatory.
%   The only trick lies in the calculation of C1 and C2. To be brief, this
%   is essentially the two-dimensional analog of the condition that needs
%   to be satisfied by a function F(x) that has a zero in the interval
%   [a,b], namely
%           F(a)*F(b) <= 0
%   C1 and C2 exactly do this for each segment of curves 1 and 2
%   respectively. If this condition is satisfied simultaneously for two
%   segments then we know that they will cross at some point. 
%   Each factor of the 'C' arrays is essentially a matrix containing 
%   the numerators of the signed distances between points of one curve
%   and line segments of the other.
    %...Argument checks and assignment of L2
    error(nargchk(1,2,nargin));
    if nargin == 1,
        L2 = L1;    hF = @lt;   %...Avoid the inclusion of common points
    else
        L2 = varargin{1}; hF = @le;
    end
       
    %...Preliminary stuff
    x1  = L1(1,:)';  x2 = L2(1,:);
    y1  = L1(2,:)';  y2 = L2(2,:);
    dx1 = diff(x1); dy1 = diff(y1);
    dx2 = diff(x2); dy2 = diff(y2);
    
    %...Determine 'signed distances'   
    S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
    S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);
    
    C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
    C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';
    %...Obtain the segments where an intersection is expected
    [i,j] = find(C1 & C2); 
    if isempty(i),P = zeros(2,0);return; end;
    
    %...Transpose and prepare for output
    i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
    L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
    i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0
    
    %...Solve system of eqs to get the common points
    P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
                dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';
              
    function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
    end
end

