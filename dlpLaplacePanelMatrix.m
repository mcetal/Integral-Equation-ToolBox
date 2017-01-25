function DLP = dlpLaplacePanelMatrix(nPanel, npt, w, z, ds, Nz, kappa)
% DLPLAPLACEPANELMATRIX(alpha, nPanel, npt, T, w, W, z, ds, Nz, kappa) 
%  Constructs system matrix for Fredholm IE for the Dirichlet BVP for 
%  Laplace: solution is represented by a DLP.
%  Panel based quadrature.
%
% INPUTS:
%   nPanel:
%       Number of panels
%   npt:
%       Number of nodes per panel 
%   t:
%       Parameter value at nodes
%   T:
%       Canonical Gauss-Legendre nodes
%   w:
%       Regular quadrature weights
%   W: 
%       Canonical quadrature weights
%   z:
%       Boundary points
%   ds:
%       |dz/dt|
%   Nz:
%       Outward pointing normal
%   kappa:
%       Curvature
%
% OUTPUTS:
%   DLP:
%       System matrix for IE

    npts = nPanel*npt;
        
%
% Construct system matrix
    
    DLP = zeros(npts, npts);
    
    for iPanel = 1: nPanel
        i = (iPanel - 1)*npt + (1: npt);
        
        Diagonal = 1 + 0.5*ds(i).*w.*kappa(i)/pi;
        
        for index = 1: npt
            itarg = (iPanel - 1)*npt + index;
            
%
%        off diagonal panels
            for jPanel = 1: nPanel
                j = (jPanel - 1)*npt + (1: npt);
                dR = z(j) - z(itarg);       
                gradU = dR./abs(dR).^2;
                DLP(itarg, j) ...
                    = w.*(real(gradU.*conj(Nz(j)))).*ds(j)/pi;
            end
%
%        diagonal 
            DLP(itarg, itarg) = Diagonal(index);
            
        end
    end
    
end

