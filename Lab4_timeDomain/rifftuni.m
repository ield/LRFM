function realsignal=rifftuni(unilateralspectrum)
% RIFFTUNI  Real Inverse of a Unilateral Fast Fourier Transform
%   realsignal = RIFFTUNI(unilateralspectrum) computes a real signal form
%                                           its unilateral spectrum (only 
%                                           positive frequencies).
%
%   See also fftuni, initializeFT

%   J. Esteban 7/10/2016

M=length(unilateralspectrum);
N=2*(M-1);
unilateralspectrum(1)=2*unilateralspectrum(1);
unilateralspectrum=(M-1)*unilateralspectrum;
realsignal=real(ifft([unilateralspectrum; ... 
                     flipud(conj(unilateralspectrum(2:N/2)))]));