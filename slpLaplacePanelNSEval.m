function [val] = slpLaplacePanelNSEval(nPanel, npt, w, z, zP,  LGammaP, ...
                                       dz, Nz, ds, rho, zTarg, AV)
% SLPLAPLACEPANELNSEVAL(alpha, nPanel, npt, w, z,zP, Nz, dz, ds, rho, zTarg) 
%  Evaluate the double layer potential with density rho and the point
%  zTarg with product integration activated (i.e. zTarg is close to
%  source). See section 6 of Helsing & Holst (2015)
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
%   AV:
%       Matrix to compute average density - used as regularization for
%       SLP
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
    [d2GammaP, iclose] = min(abs(zTarg - z));
    iPanel = floor((iclose-1)/npt) + 1; 
    
    val = slpLaplacePanelEval(nPanel, npt, w, z, ds, rho, zTarg, AV);
    
    GL = 0.5/pi;

%
% Gamma_p
    dstol = tolfac*LGammaP(iPanel);
    if d2GammaP < dstol
        j = (iPanel-1)*npt + (1:npt);
        za = zP(iPanel);
        zb = zP(iPanel+1);
        rpwj = w.*dz(j);
    
        [wcorrL, ~] = wLCinit(za, zb, zTarg, z(j), Nz(j), rpwj, npt);
    
        dval = sum(GL.*rho(j).*ds(j).*w.*wcorrL);
        val = val + dval;
    end
%
% Gamma_p-1
    jPanel = iPanel-1;
    if jPanel == 0
        jPanel = nPanel;
    end
    j = (jPanel-1)*npt + (1:npt);
    d2GammaP = min(abs(zTarg - z(j)));
    dstol = tolfac*LGammaP(jPanel);
    if d2GammaP < dstol
        za = zP(jPanel);
        zb = zP(jPanel+1);
        rpwj = w.*dz(j);
    
        [wcorrL, ~] = wLCinit(za, zb, zTarg, z(j), Nz(j), rpwj, npt);
    
        dval = sum(GL.*rho(j).*ds(j).*w.*wcorrL);
        val = val + dval;
    end
%
% Gamma_p+1
    jPanel = iPanel+1;
    if jPanel == nPanel + 1
        jPanel = 1;
    end
    j = (jPanel-1)*npt + (1:npt);
    d2GammaP = min(abs(zTarg - z(j)));
    dstol = tolfac*LGammaP(jPanel);
    if d2GammaP < dstol
        za = zP(jPanel);
        zb = zP(jPanel+1);
        rpwj = w.*dz(j);
    
        [wcorrL, ~] = wLCinit(za, zb, zTarg, z(j), Nz(j), rpwj, npt);
    
        dval = sum(GL.*rho(j).*ds(j).*w.*wcorrL);
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

