function out=gate(t,tini,tfin,tram)
% GATE  Rectangular time gate with sinusoidal edges
%   out=GATE(t,tini,tfin,tram) computes a rectangular time gate with 
%                              sinusoidal rising and falling edges. The
%                              gate begins at tini and ends at tfin, and 
%                              has unity value from tini+tram to tfin-tram.
%
%   t = time vector
%   tini = initial gate time
%   tfin = final gate time
%   tram = time duration of sinusoidal edge
%
%   out = column vector with gate values

%   J. Esteban 01/02/2018

out=1*((t>=tini+tram)&(t<=tfin-tram))+...
    (sin(pi*(t-tini)/(2*tram)).^2).*((t>tini)&(t<tini+tram))+...
    (sin(pi*(t-tfin)/(2*tram)).^2).*((t>tfin-tram)&(t<tfin));
