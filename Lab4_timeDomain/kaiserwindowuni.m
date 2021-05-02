function out=kaiserwindowuni(N,beta)
% KAISERWINDOWUNI  Unilateral Kaiser-Bessel window with parameter beta
%   out=KAISERWINDOWUNI(N,beta) computes the unilateral (only positive 
%                                           frequencies) Kaiser-Bessel 
%                                           window with parameter beta.
%
%   N = Number of frequency points
%   beta = window paramenter
%
%   out = column vector with window values

%   J. Esteban 01/02/2018

out=besseli(0,beta*sqrt(1-((0:N)'/N).^2))/besseli(0,beta);