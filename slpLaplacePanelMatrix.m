function SLP = slpLaplacePanelMatrix(nPanel, npt, t, T, w, W, ...
                                               z, ds)

% SLPLAPLACEPANELMATRIX(alpha, nPanel, npt, T, w, W, z, ds, Nz, kappa) 
%  Constructs system matrix for single layer potential.
%  Quadrature is kernel-split, panel-based a la Helsing & Holst (2015).
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
%
% OUTPUTS:
%   SLP:
%       System matrix for IE
%

    npts = nPanel*npt;
    dt = 2*pi/nPanel;
    
    if npt == 16
        ActivateProdFac = 1;
    else
        ActivateProdFac = 0.7;
    end
%
% Construct corrected quadrature weights; see equations (30) & (31) of 
% Helsing & Holst
    
    WPm1 = WfrakLinit(-2, 1, T, npt);
    WP = WfrakLinit(0, 1, T, npt);
    WPp1 = WfrakLinit(2, 1, T, npt);
    
    wCorrectDiag = bsxfun(@rdivide, WP, W');
    wCorrectSub = bsxfun(@rdivide, WPp1, W');
    wCorrectSup = bsxfun(@rdivide, WPm1, W');
    
    for i = 1: npt
        row = log(abs(T(i) - T))';
        row(i) = 0;
        wCorrectDiag(i, :) = wCorrectDiag(i, :) - row;
        row = log(abs(T(i) + 2 - T))';
        wCorrectSub(i, :) = wCorrectSub(i, :) - row;
        row = log(abs(T(i) - 2 - T))';
        wCorrectSup(i, :) = wCorrectSup(i, :) - row;
    end
    GL = 0.5/pi;
    
%
% Construct system matrix
    M0 = zeros(npts, npts);
    MS = zeros(npts, npts);
    
    for iPanel = 1: nPanel
        iPm1 = mod(iPanel - 2, nPanel) + 1;
        iPp1 = mod(iPanel, nPanel) + 1;
        i = (iPanel - 1)*npt + (1: npt);
        
        
        for index = 1: npt
            itarg = (iPanel - 1)*npt + index;
            
%
% regular quadrature for distant panels
            for jp = 1: nPanel - 3
                jPanel = mod(iPp1 + jp -1, nPanel) + 1;
%            disp(['    jPanel = ', num2str(jPanel)])
                j = (jPanel - 1)*npt + (1: npt);
                M0(itarg, j) = (Kernel(z(j), z(itarg)).*ds(j).*w)';
            end
%
% Construct M^*
%    Diagonal Block
            jPanel = iPanel;
%        regular quadrature
            j = (jPanel - 1)*npt + (1: npt);
            M0(itarg, j) = (Kernel(z(j), z(itarg)).*ds(j).*w)';
            M0(itarg, itarg) = 0;
%            M0(itarg, itarg) = 0.5*ds(itarg)*w(index)/pi;
            
%       correction            
            wCorrect = wCorrectDiag(index, :)';
            wCorrect(index) = wCorrectDiag(index, index) ...
                                  + log(0.5*ds(i(index))*dt);
            MS(itarg, j) = (GL.*ds(j).*w.*wCorrect)';
            
%       
%    Sub Diagonal Block
            jPanel = iPm1;
            ta = (jPanel-1)*dt;
            tb = ta + dt;
            tmid = 0.5*(ta + tb);
            if jPanel == nPanel
                tmid = tmid - 2*pi;
            end
            
%        regular quadrature
            j = (jPanel - 1)*npt + (1: npt);
            M0(itarg, j) = (Kernel(z(j), z(itarg)).*ds(j).*w)';
            
%        correction            
            if abs(t(itarg) - tmid) < ActivateProdFac*dt
                MS(itarg, j) = (GL.*ds(j).*w.*wCorrectSub(index,:)')';
            end
%    Super Diagonal Block
            jPanel = iPp1;
            ta = (jPanel-1)*dt;
            tb = ta + dt;
            tmid = 0.5*(ta + tb);
            if jPanel == 1
                tmid = tmid + 2*pi;
            end
            
%        regular quadrature
            j = (jPanel - 1)*npt + (1: npt);
            M0(itarg, j) = (Kernel(z(j), z(itarg)).*ds(j).*w)';
            
%        correction            
            if abs(t(itarg) - tmid) < ActivateProdFac*dt
                MS(itarg, j) = (GL.*ds(j).*w.*wCorrectSup(index,:)')';
            end
        end
        
        
    end
    
    SLP = M0 + MS;

end

function [Kern] = Kernel(zSource, zTarget)
    dR = zSource - zTarget;
    Kern = 0.5*(log(abs(dR))) /pi;       
end

function WfrakL = WfrakLinit(trans, scale, tfrak, npt)
%
% See Appendix A, Helsing & Holst

    A = fliplr(vander(tfrak));
    tt = trans + scale*tfrak;
    Q = zeros(npt);
    p = zeros(1, npt+1);
    c=(1-(-1).^(1:npt))./(1:npt);
    for m = 1: npt
        p(1) = log(abs((1-tt(m))/(1+tt(m))));
        p1 = log(abs(1-tt(m)^2));
        for k = 1:npt
            p(k+1)=tt(m)*p(k)+c(k);
        end
        Q(m, 1: 2: npt-1) = p1 - p(2: 2: npt);
        Q(m, 2: 2: npt) = p(1) - p(3: 2: npt+1);
        Q(m, :) = Q(m, :) ./ (1:npt);
    end
    
    WfrakL=Q/A;
end
  