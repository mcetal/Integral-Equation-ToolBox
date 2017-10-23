function [val] = dlpLaplacePanelNSEval(nPanel, npt, nBody, w, z, zP,  ...
                                      LGammaP, Nz, dz, ds, rho, zTarg, ...
                                      iTarg)
% DLPYUKAWAPANELNSEVAL(alpha, nPanel, npt, w, z,zP, Nz, dz, ds, rho, zTarg) 
%  Evaluate the double layer potential with density rho and the point
%  zTarg with product integration activated (i.e. zTarg is close to
%  source). See section 6 of Helsing & Holst (2015)
%
% INPUTS:
%   alpha:
%       Parameter in PDE: u - \alpha^2 \Delta u = 0
%   nPanel:
%       Number of panels
%   npt:
%       Number of nodes per panel 
%   nBody:
%       Number of component curves
%   w:
%       Regular quadrature weights
%   z:
%       Boundary points
%   zP:
%       End points of panels on curve
%   LGammaP:
%       arc length of each panel
%   Nz:
%       Outward pointing normal
%   dz:
%       Derivative of z wrt t
%   ds:
%       |dz/dt|
%   rho:
%       Double layer density
%   zTarg:
%       Target point
%   iTarg:
%       Component curve zTarg is close to
%
% OUTPUTS:
%   val:
%       value of DLP

    if npt ==16
        tolfac = 1.1;
    else
        tolfac = 0.3;
    end

%
% Find closest panel
    [d2GammaP, iclose] = min(abs(zTarg - z(:, iTarg)));
    iPanel = floor((iclose-1)/npt) + 1; 
    
    val = dlpLaplacePanelEval(nPanel, npt, nBody, w, z, Nz, ds, rho, ...
                              zTarg);
    
%
% Gamma_p
    dstol = tolfac*LGammaP(iPanel, iTarg);
    if d2GammaP < dstol
        j = (iPanel-1)*npt + (1:npt);
        za = zP(iPanel, iTarg);
        zb = zP(iPanel+1, iTarg);
        rpwj = w.*dz(j, iTarg);
    
        [~, wcmpC] = wLCinit(za, zb, zTarg, z(j, iTarg), Nz(j, iTarg), ...
                             rpwj, npt);
    
        dval = sum(rho(j, iTarg).*wcmpC)/(2*pi);
        val = val + dval;
    end
%
% Gamma_p-1
    jPanel = iPanel-1;
    if jPanel == 0
        jPanel = nPanel;
    end
    j = (jPanel-1)*npt + (1:npt);
    d2GammaP = min(abs(zTarg - z(j, iTarg)));
    dstol = tolfac*LGammaP(jPanel, iTarg);
    if d2GammaP < dstol
        za = zP(jPanel, iTarg);
        zb = zP(jPanel+1, iTarg);
        rpwj = w.*dz(j, iTarg);
    
        [~, wcmpC] = wLCinit(za, zb, zTarg, z(j, iTarg), Nz(j, iTarg), ...
                             rpwj, npt);
    
        dval = sum(rho(j, iTarg).*wcmpC)/(2*pi);
        val = val + dval;
    end
%
% Gamma_p+1
    jPanel = iPanel+1;
    if jPanel == nPanel + 1
        jPanel = 1;
    end
    j = (jPanel-1)*npt + (1:npt);
    d2GammaP = min(abs(zTarg - z(j, iTarg)));
    dstol = tolfac*LGammaP(jPanel, iTarg);
    if d2GammaP < dstol
        za = zP(jPanel, iTarg);
        zb = zP(jPanel+1, iTarg);
        rpwj = w.*dz(j, iTarg);
    
        [~, wcmpC] = wLCinit(za, zb, zTarg, z(j, iTarg), Nz(j, iTarg), ...
                             rpwj, npt);
    
        dval = sum(rho(j, iTarg).*wcmpC)/(2*pi);
        val = val + dval;
    end
    
end

function [wcorrL,wcmpC, A] = wLCinit(ra, rb, r, rj, nuj, rpwj, npt)
% See appendix B, Helsing & Holst
% *** ztgtr is target vector in transformed plane ***
    dr = (rb-ra)/2;
    rtr = (r - (rb+ra)/2)/dr;
    rjtr = (rj - (rb+ra)/2)/dr;
    A = fliplr(vander(rjtr)).';
    p = zeros(npt+1, 1);
    q = zeros(npt, 1);
    c=(1-(-1).^(1:npt))./(1:npt);
    p(1) = log(1 - rtr) - log(-1 - rtr);
    p1 = log(1 - rtr) + log(-1 - rtr);
    if imag(rtr) < 0 && abs(real(rtr)) < 1
        p(1) = p(1) + 2*1i*pi;
        p1 = p1 - 2*1i*pi;
    end
    for k = 1:npt
        p(k+1) = rtr*p(k)+c(k);
    end
    q(1:2:npt-1) = p1 - p(2:2:npt);
    q(2:2:npt) = p(1) - p(3:2:npt+1);
    q = q./(1:npt)';
    wcorrL = imag(A\q*dr.*conj(nuj))./abs(rpwj)-log(abs((rj-r)/dr));  
    wcmpC = imag(A\p(1:npt)-rpwj./(rj-r));
    
    
 end

