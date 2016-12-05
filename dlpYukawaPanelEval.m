function [val] = dlpYukawaPanelEval(alpha, nPanel, npt, w, z, Nz, ds, ...
                                    sigma, zTarg)
% DLPYUKAWAPANELEVAL(alpha, nPanel, npt, w, z, Nz, ds, sigma, zTarg) 
%  Evaluate the double layer potential with density sigma and the point
%  zTarg using composite Gauss-Legendre.
%
% INPUTS:
%   alpha:
%       Parameter in PDE: u - \alpha^2 \Delta u = 0
%   nPanel:
%       Number of panels
%   npt:
%       Number of nodes per panel 
%   w:
%       Regular quadrature weights
%   z:
%       Boundary points
%   Nz:
%       Outward pointing normal
%   ds:
%       |dz/dt|
%   sigma:
%       Double layer density
%   zTarg:
%       Target point
%
% OUTPUTS:
%   val:
%       value of DLP

    val = 0;
    for jpanel = 1: nPanel
        j = (jpanel-1)*npt + (1: npt);
        Kern = Kernel(alpha, z(j), zTarg, Nz(j));
        val = val - sum(Kern.*ds(j).*sigma(j).*w)/(2*alpha^2);
    end
end 

function [Kern] = Kernel(alpha, zSource, zTarget, Nz)
    dR = zSource - zTarget;
    gradU = besselk(1, abs(dR)/alpha).*dR./abs(alpha*dR);
    Kern = real(gradU.*conj(Nz)) /pi;       
end
