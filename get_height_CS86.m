% GET_HEIGHT_CS86
%
%   [height, wind_speed] = GET_HEIGHT_CS86(dw_range, cw_range, diameter, density)
%
%   Calculates plume height (km above sampling height) and wind speed (m s-1)
%   using an implementation of the model of Carey and Sparks (1986) based on
%   polynomial fits. 4 inputs are required:
%       - Isopleth downwind range (km)
%       - Isopleth crosswind range (km)
%       - Clast diameter (cm)
%       - Clast density (kg m-3)
%
%   Pass the arguments in the command line to use this script as a function. 
%   Passing no argument will start the GUI.
%
function [height, wind] = get_height_CS86(varargin)
% This function is part of the TError package by S. Biass (Dept of Earth
% Sciences, University of Geneva), G. Bagheri, W. Aeberhard and C.
% Bonadonna.
% - http://dx.doi.org/10.5038/2163-338X.1.2 
% - https://vhub.org/resources/3701
%
% HISTORY
%   2015/10/04: First release
%
% This script was written by S. Biass and G. Bagheri
% Copyright: S. Biass and G. Bagheri (University of Geneva) on the script
% This script is distributed under a GNU GPL3 license.
% For any question: sebastien.biasse@unige.ch

%{
This is a free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    It is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with it. If not, see <http://www.gnu.org/licenses/>.
%}

% Process input parameters
if nargin == 0          % When no input is specified, a GUI opens
    prompt      = {'Downwind range (km):','Crosswind range (km):', 'Clast diameter (cm):', 'Clast density (kg m-3):'};
    dlg_title   = 'Carey & Sparks 1986';
    num_lines   = 1;
    answer      = inputdlg(prompt,dlg_title,num_lines);
    if isempty(answer)
        return
    end
    % Retrieve inputs:
    dw          = str2double(answer(1)); % Crosswind range
    cw          = str2double(answer(2)); % Downwind range
    d           = str2double(answer(3)); % Diameter
    den         = str2double(answer(4)); % Density
elseif nargin >0 && nargin ~= 4
    error('myapp:argChk', 'Wrong number of arguments. Use either 0 arguments to start the GUI or 4 arguments to use the function. Type\n>> help get_height_CS86\nfor more help.')   
else
    dw          = varargin{1}; % Crosswind range
    cw          = varargin{2}; % Downwind range
    d           = varargin{3}; % Diameter
    den         = varargin{4}; % Density
end

d = d/100; % Convert the clast diameter to metres

% Setup iteration on elongation and flatness
fl_min = 0.1;       % Minimum flatness
fl_max = 1;         % Maximum flatness
el_min = 0.1;       % Minimum elongation
el_max = 1;         % Maximum elongation

shape_inc = 0.01;   % Increments
fl_vec = fl_min:shape_inc:fl_max;
el_vec = el_min:shape_inc:el_max;

% Storage matrices
stor_ht = zeros(length(el_vec), length(fl_vec));
stor_wd = zeros(length(el_vec), length(fl_vec));
stor_deq=zeros(length(el_vec), length(fl_vec));

% Iteration on elongation and flatness
for iF = 1:length(fl_vec)
    for iE = 1:length(el_vec)
        vtc             = get_velocity(d, den, fl_vec(iF), el_vec(iE));     % Get the terminal velocity of a given flatness/elongation
        d_eq            = (vtc/1.054)^2 * 1.096 / (9.81*den);               % Matching terminal velocity with equation (3) of Carey and Sparks
        [Wtmp, Htmp]    = interp_height(dw, cw, d_eq, den);                 % Retrieve height and wind for d_eq
        % Fill storage
        stor_ht(iE,iF)  = Htmp;
        stor_wd(iE,iF)  = Wtmp;
        stor_deq(iE,iF) = d_eq;
    end
end


% For each elongation, finds flatness that gives minimum difference of
% terminal velocity (ie. diameter) between Carey and Sparks (1986) and drag
% coefficient of Bagheri and Bonadonna (submitted).
plotX = zeros(length(fl_vec));
plotY = zeros(length(el_vec));

for iY = 1:length(el_vec)
    [~, iX] = min(abs(stor_deq(iY,:)-d));
    plotX(iY)   = fl_vec(iX)+0.5*shape_inc;
    plotY(iY)   = el_vec(iY)+0.5*shape_inc;
end


% Plot plume height
figure; 
[C,h] = contourf(fl_vec, el_vec, stor_ht); xlabel('Flatness'); ylabel('Elongation'); axis equal; hold on;
axis([.1 1 .1 1]);
rectangle('Position', [0.4,0.4,1,1],'LineWidth', 2);
title(['Plume height - ', num2str(d*100, '%1.1f'), ' cm, ', num2str(den), ' kg m^-^3']);
clabel(C,h,'LabelSpacing', 1000);
cb = colorbar;
ylabel(cb, 'Plume height (km)');

plot(plotX, plotY, 'xr')

% Plot wind speed
figure; 
[C,h] = contourf(fl_vec, el_vec, stor_wd); xlabel('Flatness'); ylabel('Elongation'); axis equal; hold on;
axis([.1 1 .1 1]);
rectangle('Position', [0.4,0.4,1,1],'LineWidth', 2);
title(['Wind speed - ', num2str(d*100, '%1.1f'), ' cm, ', num2str(den), ' kg m^-^3']);
clabel(C,h,'LabelSpacing', 1000);
cb = colorbar;
ylabel(cb,'Wind ppeed (m s-1)');

plot(plotX, plotY, 'xr')



% Display final results
[wind, height] = interp_height(dw, cw, d, den);
display(sprintf('Estimates from Carey and Sparks (1986):\n - Plume height:\t %.1f km above sampling height\n - Wind speed:\t\t %.0f m/s\n', height, wind))

% Get the terminal velocity of a given flatness/elongation (Bagheri and
% Bonadonna, submitted)
function vtc = get_velocity(d, den, el, fl)
mass    = den*pi*d^3/6;

vtc     = 10^-12;                       % Initial terminal velocity
acc     = 9.81;                         % Particle acceleration
den_f   = 1.096;                        % Fluid density
vis_f   = 1.8975*10^-5;                 % Fluid dynamic viscosity
Fs      = fl*el^1.3;                    % Stokes shape descriptor
ks      = 0.5*(Fs^(1/3) + Fs^(-1/3));   % Stokes drag correction
Fn      = fl^2*el;                      % Newton shape descriptor
kn      = 10^(.45* (-log10(Fn))^.99);   % Newtown drag correction
dT      = 0.05;                         % Time step for falling simulation

while acc > 0.01*9.81                   % Criteria for reaching terminal velocity
    Re      = den_f*vtc*d/vis_f;        % Reynold's number
    Re_star = Re*kn/ks;                 % Normalized Reynold's number
    CD      = kn * (24/Re_star * (1 + 0.125 * Re_star^(2/3)) + 0.46/ (1+5330/Re_star)); % Drag coefficient
    F_drag  = 0.5*den_f*0.25*pi*d^2*CD*vtc^2; % Drag force
    acc     = (mass*9.81-F_drag)/mass;  % Update acceleration
    vtc     = acc*dT + vtc;             % Update falling velocity
end


function [wind, height] = interp_height(dw, cw, d, den)



%% 1. Find values of height and wind for a given set of cw/dw on all Figures 16 of Carey and Sparks (1986)
% Figure 16a
d_a=0.8e-2;
den_a=2500;

a0=dw;
% Polynomial fit
a10=-9e-5*dw^3+0.0117*dw^2+.4887*dw;
a20=-8e-5*dw^3+0.0129*dw^2+.252*dw;
a30=-6e-5*dw^3+0.0114*dw^2+.1229*dw;

if (cw>=a0)
    wind_a=0;
elseif (cw<a0) && (cw>=a10)
    wind_a=0.+(10-0)*(cw-a0)/(a10-a0);
elseif (cw<a10) && (cw>=a20)
    wind_a=10.+(20-10)*(cw-a10)/(a20-a10);
elseif (cw<a20)
    wind_a=20.+(30-20)*(cw-a20)/(a30-a20); 
end
height_a=(-0.01977*cw^3 + 1.524*cw^2 + 16.26*cw + 0.06977) / (cw + 1.928);

% Figure 16b
d_b=1.6e-2;
den_b=2500;
b0=dw;
% Polynomial fit
b10=-2e-5*dw^3+0.0092*dw^2+.5168*dw;
b20=-10e-5*dw^3+0.0151*dw^2+.2577*dw;
b30=-10e-5*dw^3+0.0153*dw^2+.1191*dw;

if (cw>=b0)
    wind_b=0;
elseif (cw<b0) && (cw>=b10)
    wind_b=0.+(10-0)*(cw-b0)/(b10-b0);
elseif (cw<b10) && (cw>=b20)
    wind_b=10.+(20-10)*(cw-b10)/(b20-b10);
elseif (cw<b20)
    wind_b=20.+(30-20)*(cw-b20)/(b30-b20); 
end
height_b=(-0.02357*cw^3 + 1.672*cw^2 + 17.69*cw + 0.08577) / (cw + 1.557);   

% Figure 16c
d_c=3.2e-2;
den_c=2500;
c0=dw;
% Rational fit
c10=(0.01254 *dw^5 + 0.2165*dw^4 -9.467 *dw^3 +  85.12*dw^2 -199.4*dw -15.95) /(dw^3 -23.58*dw^2 +  183.3*dw -415);
c20= (0.01259*dw^5 +  0.0369 *dw^4 -5.405*dw^3 +55.24*dw^2 -104.7*dw + 22.13) /(dw^3 -23.97*dw^2 + 183.7*dw -294.2);
c30= (-0.0314*dw^5 + 13.41*dw^4  -85.39*dw^3  -4602*dw^2 + 5.544e+04*dw -3507 ) /(dw^3 + 1061*dw^2  -3.068e+04*dw + 2.413e+05);

if (cw>=c0)
    wind_c=0;
elseif (cw<c0) && (cw>=c10)
    wind_c=0.+(10-0)*(cw-c0)/(c10-c0);
elseif (cw<c10) && (cw>=c20)
    wind_c=10.+(20-10)*(cw-c10)/(c20-c10);
elseif (cw<c20)
    wind_c=20.+(30-20)*(cw-c20)/(c30-c20); 
end
height_c=(-0.03649*cw^3 +2.188*cw^2 + 15.27*cw -0.03349) / (cw +   0.8307);

% Figure 16d
d_d=6.4e-2;
den_d=2500;
d0=dw;
% Polynomial fit
d10=-0.0009*dw^3+0.0379*dw^2+.4446*dw;
d20=-0.0007*dw^3+0.04*dw^2+.1919*dw;
d30=-0.0006*dw^3+0.0354*dw^2+.0867*dw;

if (cw>=d0)
    wind_d=0;
elseif (cw<d0) && (cw>=d10)
    wind_d=0.+(10-0)*(cw-d0)/(d10-d0);
elseif (cw<d10) && (cw>=d20)
    wind_d=10.+(20-10)*(cw-d10)/(d20-d10);
elseif (cw<d20)
    wind_d=20.+(30-20)*(cw-d20)/(d30-d20); 
end
height_d=(-0.05547*cw^3 + 2.729*cw^2 + 18.71*cw + 0.001767) / (cw + 0.9346);


%% 2. Interpolates the values of height and wind for a given diameter/density
if (d*den<=d_a*den_a)
    height=height_a+(height_b-height_a)*(d*den-d_a*den_a)/(d_b*den_b-d_a*den_a);
    wind=wind_a+(wind_b-wind_a)*(d*den-d_a*den_a)/(d_b*den_b-d_a*den_a);
elseif (d*den>=d_a*den_a)&& (d*den<d_b*den_b)
    height=height_a+(height_b-height_a)*(d*den-d_a*den_a)/(d_b*den_b-d_a*den_a);
    wind=wind_a+(wind_b-wind_a)*(d*den-d_a*den_a)/(d_b*den_b-d_a*den_a);
elseif (d*den>=d_b*den_b)&& (d*den<d_c*den_c)
    height=height_b+(height_c-height_b)*(d*den-d_b*den_b)/(d_c*den_c-d_b*den_b);
    wind=wind_b+(wind_c-wind_b)*(d*den-d_b*den_b)/(d_c*den_c-d_b*den_b);
elseif (d*den>=d_c*den_c)   
    height=height_c+(height_d-height_c)*(d*den-d_c*den_c)/(d_d*den_d-d_c*den_c);
    wind=wind_c+(wind_d-wind_c)*(d*den-d_c*den_c)/(d_d*den_d-d_c*den_c);
end
