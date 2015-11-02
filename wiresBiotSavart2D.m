% 2D Biot Savart simulation of wires. Two parallel wires in x direction. Separated by
% 1.4mm in the y direction. Carrying opposite currents. Plots the field
% in the y-z plane.
%Debi Prasad Pattnaik; U of Nottingham
clear all;
close all;
clc;

mu = pi*4e-7;
I = 10;
L = 20e-3;

NumCurrentElements = 1000;
NumPoints = 300;


% Define the grid.
x=linspace(-2e-3,2e-3,NumPoints);
y=linspace(-1e-3,1e-3,NumPoints);
z=linspace(0.1e-3,4e-3,NumPoints);
plane = 'yz';
if strcmp('xy',plane)
    [xgrid,ygrid] = meshgrid(x,y);
    zgrid = 0;
elseif strcmp('xz',plane)
    [xgrid,zgrid] = meshgrid(x,z);
    ygrid = 0;
elseif strcmp('yz',plane)
    [ygrid,zgrid] = meshgrid(y,z);
    xgrid = 0;
end


% Defines the co-ordinates of the lines of charge.
zstepsize(1) = L/NumCurrentElements;
v(1,:) = [1 0 0];
xsource(1,:) = linspace(-L/2,L/2,NumCurrentElements);
ysource(1,:) = zeros(1,NumCurrentElements)-0.7e-3;
zsource(1,:) = zeros(1,NumCurrentElements);

zstepsize(2) = L/NumCurrentElements;
v(2,:) = [-1 0 0];
xsource(2,:) = linspace(-L/2,L/2,NumCurrentElements);
ysource(2,:) = zeros(1,NumCurrentElements)+0.7e-3;
zsource(2,:) = zeros(1,NumCurrentElements);


% Finds the magnetic field contribution for each element of current along
% the rings.

Bxtot = 0; Bytot = 0; Bztot = 0;
for wire = 1:2

    for n = 1:NumCurrentElements
        
        rx = xgrid - xsource(wire,n);
        ry = ygrid - ysource(wire,n);
        rz = zgrid - zsource(wire,n);
        r = sqrt(rx.^2 + ry.^2+rz.^2);
        
        dI = zstepsize(wire) * v(wire,:);
        xcrossed = (dI(2)*rz) - (ry*dI(3));
        ycrossed = (dI(3)*rx) - (rz*dI(1));
        zcrossed = (dI(1)*ry) - (rx*dI(2));
        
        Bx = (mu/(4*pi))*I*xcrossed./(r.^3);
        By = (mu/(4*pi))*I*ycrossed./(r.^3);
        Bz = (mu/(4*pi))*I*zcrossed./(r.^3);
        
        Bxtot = Bxtot+Bx;
        Bytot = Bytot+By;
        Bztot = Bztot+Bz;
        
        
    end
wire
end


Bmag = sqrt(Bxtot.^2 + Bytot.^2 + Bztot.^2);
figure(1)
contour(y,z,log(Bmag),200);



 figure(2)
 iplotline = findindex(y,0);
 plot(z,Bmag(:,iplotline));
 figure(3)
 dB = gradient(Bmag(:,iplotline),z);
 plot(z,dB*10000/100)
title('Static Magnetic Field Gradient G/cm');
 xlabel('z-direction');
 ylabel('dB');



% Btheory = 2e-7 * I./(z);
% figure(3)
% plot(z,Btheory)