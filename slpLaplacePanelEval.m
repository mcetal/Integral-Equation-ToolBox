function [val] = slpLaplacePanelEval(nPanel, npt, w, z, ds, sigma, ...
                                     zTarg)
% DLPLAPLACEPANELEVAL(nPanel, npt, w, z, Nz, ds, sigma, zTarg) 
%  Evaluate the double layer potential with density sigma and the point
%  zTarg using composite Gauss-Legendre.
%
% INPUTS:
%   nPanel:
%       Number of panels
%   npt:
%       Number of nodes per panel 
%   w:
%       Regular quadrature weights
%   z:
%       Boundary points
%   ds:
%       |dz/dt|
%   sigma:
%       Double layer density
%   zTarg:
%       Target point
%
%
% OUTPUTS:
%   val:
%       value of DLP

    val = 0;
    for jpanel = 1: nPanel
        j = (jpanel-1)*npt + (1: npt);
        dR = z(j) - zTarg;       
%        kernel = 0.5*(log(abs(dR)) + 1).*ds(j)/pi; 
        kernel = 0.5*(log(abs(dR)) ).*ds(j)/pi; 
        val = val + sum(w.*kernel.*sigma(j));
    end
end

