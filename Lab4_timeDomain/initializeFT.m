function [t,f]=initializeFT(n_frec,Delta_f)
% INITIALIZEFT defines time and frequency column vectors
%   [t,f]=INITIALIZEFT(n_frec,Delta_f)
%
%   n_frec = Number of frequencies, increased up to a power of 2
%   Delta_f = Frequency increment
%
%   t = Time column vector:  [0; Delta_t; 2Delta_t; 3Delta_t;... t_MAX]
%   f = Frequency column vector: [0; Delta_f; 2Delta_f; 3Delta_f;... f_MAX]
%
%   See also fftuni, rifftuni

%   J. Esteban 7/10/2016

N=1+ceil(log(n_frec)/log(2));
M=2^N;

f_MAX=Delta_f*M/2;
f=linspace(0,f_MAX,M/2+1)';

Delta_t=1/(2*f_MAX);
t_MAX=Delta_t*(M-1);
t=linspace(0,t_MAX,M)';
