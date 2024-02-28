
% The Matlab/Octave file used to generate the image in Fig. 1 of Chapter 19
% (Keener and Sneyd, Mathematical Physiology)

clear all
close all 
clc

mach = zeros(1000,2000);  % A nice big matrix, to get decent resolution

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
 