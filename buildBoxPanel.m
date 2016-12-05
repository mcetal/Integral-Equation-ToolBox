function [xBox, yBox, igrid, LGammaP] = ...
                    buildBoxPanel(M, nPanel, npt, w, z, dz, ds)
% BUILDBOXPANEL(M, nPanel, npt, w, z, dz, ds)
% Embeds domain in a uniform MxM grid. Determines regular grid points - 
% i.e. points in the domain where regular quadrature is fine, and near
% singular grid points - points near the boundary where specialized
% quadrature is needed.
%
% INPUTS:
%   M:
%       Number of grid points in one dimension
%   nPanel:
%       number of panels
%   npt:
%       number of nodes per panel 
%   w:
%       Gauss Legendre Weights
%   z:
%       Boundary points
%   dz:
%       Derivative of z wrt parameter t
%   ds:
%       |dz/dt|
%
% OUTPUTS:
%   xBox, yBox:
%       Coordinates of regular grid points; MxM arrays
%   igrid:
%       MxM arrary which flags each grid point according to:
%           = 0 if outside domain
%           = 1 if regular grid point inside domain
%           = 2 if near-singular grid point
%           = -2 - temporary flag for near-singular points outside the 
%             domain
%   LGammaP:
%       arc length of each panel
%
% Determin tolerance for near-singular region according to section 6.3
% of Helsing and Holst (2015) = 1.1 |\gamma_p|

    if npt == 16
        tolfac = 1.1;
    else
        tolfac = 0.3;
    end
    
    dsMax = 0;
    arcL = 0;
    LGammaP = zeros(nPanel, 1);
    for iPanel = 1: nPanel
        i = (iPanel-1)*npt + (1: npt);
        LGammaP(iPanel) = sum(ds(i).*w);
    end
    
    dsTol = tolfac*max(LGammaP);
    
%%
% Embed domain in regular grid
    igrid = zeros(M, M);
    x0 = min(real(z));
    y0 = min(imag(z));
    xMax = max(real(z));
    yMax = max(imag(z));
    LBox = max(xMax - x0, yMax - y0);
    
    [xBox, yBox] = meshgrid(linspace(x0, x0+LBox, M), ...
                            linspace(y0, y0+LBox, M));
    z0 = x0 + 1i*y0;
    zBox = xBox + 1i*yBox;

%%
% Bin sort grid and flag points that are close
    MBin = floor(LBox/dsTol);
%
%  shift and scale so all points on box [0 MBin] x [0 MBin]
    zBoxBin = floor((zBox - z0) * MBin/LBox);
    zGammaBin = floor((z - z0) * MBin/LBox);
    
%
%  index to all boxes containing boundary points
    GammaBinIndex = unique(zGammaBin);  

% 
%  Search all boxes containing boundary points
    for iBin = 1: length(GammaBinIndex)
        ij = GammaBinIndex(iBin);
        i = real(ij);
        j = imag(ij);
        zGamma = z(zGammaBin == ij);
        dzGamma = dz(zGammaBin ==ij);
        
%
%  Flag points in neighbouring bins that are too close
        for iBoxBin = (i-1) : (i+1)
            for jBoxBin = (j-1) : (j+1)
                iCheck = find(zBoxBin == iBoxBin + 1i*jBoxBin);
                for jPoint = iCheck';
                    [dis2Gamma, iGamma] = min(abs(zBox(jPoint) - zGamma));
                    if dis2Gamma < dsTol
                        vGamma = [real(dzGamma(iGamma)), ...
                                  imag(dzGamma(iGamma)), 0];
                        zP2Gamma = zGamma(iGamma) - zBox(jPoint);
                        vP2Gamma = [real(zP2Gamma), imag(zP2Gamma), 0];
                        orientCheckVec = cross(vP2Gamma, vGamma);
                        if orientCheckVec(3) > 0 && igrid(jPoint) ~= -2
                            igrid(jPoint) = 2;
                        else
                            igrid(jPoint) = -2;
                        end
                    end
                end
             end
        end
    end

%%
% For all points not flagged, check to see if inside or outside according
% to the value of the Cauchy integral with density 1.
    tol = 1.d-12;
    iUnflagged = find(igrid==0);
    for index = iUnflagged'
        zPoint = zBox(index);
        cauchy = 0;
        for iPanel = 1: nPanel
            i = (iPanel-1)*npt + (1: npt);
            integrand = dz(i)./(z(i) - zPoint);
            cauchy = cauchy + sum(integrand.*w)/(2*pi*1i);
        end
        if abs(cauchy - 1) < tol
            igrid(index) = 1;
        end
    end
    
    igrid(igrid==-2) = 0;
    
end

    