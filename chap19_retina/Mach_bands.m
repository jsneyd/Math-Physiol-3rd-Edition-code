
%    -------------------------------------------------------------------
%
%     Mach bands.
%
%     For Chapter 19, Fig. 19.2 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%     Written by James Keener and James Sneyd.
%
%    -------------------------------------------------------------------

clear all
close all 
clc

mach = zeros(1000,2000);  % A big matrix, to get decent resolution

width = 100;    % The width of the ramp
light = 55000;   % The brightness on the bright side of the ramp
dark = 15500;     % The brightness on the dark side of the ramp

mach(:,1:1000-width)= dark;
mach(:,1000+width:2000) = light;
for i=1:2*width
    mach(:,1000-width+i) = dark + i*(light-dark)/(2*width);   % probably not the most efficient way, but good enough
end

mach = uint16(mach);
imshow(mach)
 