function DSLP = DslpDnLaplacePanelMatrix(nPanel, npt, w, z, ds, Nz, ...
                                         kappa, isingular)
% DSLPDNLAPLACEPANELMATRIX(alpha, nPanel, npt, T, w, W, z, ds, Nz, kappa) 
%  Constructs system matrix for Fredholm IE for the Neumann BVP for 
%  Laplace: solution is represented by a SLP.
%  This is the adjoint operator to the DLP.
%  Panel based quadrature.
%  System in the form [-I/2 + K/2] = f
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
%   isingular:
%       true if singular form of matrix is wanted, false if not
%
% OUTPUTS:
%   DLP:
%       System matrix for IE

    npts = nPanel*npt;
        
%
% Construct system matrix
    
    DSLP = zeros(npts, npts);
    
    for iPanel = 1: nPanel
        i = (iPanel - 1)*npt + (1: npt);
        
        Diagonal = 1 - 0.5*ds(i).*w.*kappa(i)/pi;
        
        for index = 1: npt
            itarg = (iPanel - 1)*npt + index;
            
%
%        off diagonal panels
            for jPanel = 1: nPanel
                j = (jPanel - 1)*npt + (1: npt);
                dR = z(itarg) - z(j);       
                gradU = dR./abs(dR).^2;
                DSLP(itarg, j) ...
                    = -w.*(real(gradU.*conj(Nz(itarg)))).*ds(j)/pi;
                if ~isingular
                    DSLP(itarg, j) = DSLP(itarg, j) + (w.*ds(j))'/pi;
                end
            end
%
%        diagonal 
            DSLP(itarg, itarg) = Diagonal(index);
            if ~isingular
                DSLP(itarg, itarg) ...
                    = Diagonal(index) + w(index)*ds(itarg)/pi;
            end
        end
    end
%    DSLP = -0.5*DSLP;
    
end

