function unilateralspectrum=fftuni(realsignal)
% FFTUNI  Unilateral Fast Fourier Transform
%   unilateralspectrum = FFTUNI(realsignal) computes de unilateral
%                                           spectrum (only positive 
%                                           frequencies) of a real signal.
%
%   See also rifftuni, initializeFT

%   J. Esteban 7/10/2016

N=length(realsignal);
M=N/2+1;
unilateralspectrum=2*fft(realsignal)/N;
unilateralspectrum=unilateralspectrum(1:M);
unilateralspectrum(1)=unilateralspectrum(1)/2;
