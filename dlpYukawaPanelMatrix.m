function [DLP, M0, MS] = dlpYukawaPanelMatrix(alpha, nPanel, npt, nBody,...
                                              t, T, w, W, z, ds, Nz, kappa)
% DLPYUKAWAPANELMATRIX(alpha, nPanel, npt, T, w, W, z, ds, Nz, kappa) 
%  Constructs system matrix for Fredholm IE for the Dirichlet Modified 
%  Helmholtz BVP: solution is represented by a DLP.
%  Quadrature is kernel-split, panel-based a la Helsing & Holst (2015).
%
% INPUTS:
%   alpha:
%       Parameter in PDE: u - \alpha^2 \Delta u = 0
%   nPanel:
%       Number of panels
%   npt:
%       Number of nodes per panel 
%   nBody:
%       Number of component curves in domain
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
%   M0, MS:
%       Kernel splitting of system matrices - see section 4 of 
%       Helsing & Holst (2015)

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
    
%
% Construct system matrix
    M0 = zeros(npts*nBody, npts*nBody);
    MS = zeros(npts*nBody, npts*nBody);
    
    for iBody = 1: nBody
        for iPanel = 1: nPanel
            i = (iPanel - 1)*npt + (1: npt);
            for jBody = 1: nBody
                if iBody == jBody
                    G0Diag = 0.5*ds(i, iBody).*w.*kappa(i, iBody)/pi;
                    iPm1 = mod(iPanel - 2, nPanel) + 1;
                    iPp1 = mod(iPanel, nPanel) + 1;
                    for index = 1: npt
                        itarg = (iPanel - 1)*npt + index;
            
%
% regular quadrature for distant panels
                        for jp = 1: nPanel - 3
                            jPanel = mod(iPp1 + jp -1, nPanel) + 1;
%            disp(['    jPanel = ', num2str(jPanel)])
                            j = (jPanel - 1)*npt + (1: npt);
                            M0((iBody-1)*npts + itarg, ...
                               (jBody-1)*npts + j) ...
                              = Kernel(alpha, z(j, jBody), ...
                                       z(itarg, jBody), Nz(j, jBody))...
                                       .*ds(j, jBody).*w;
                        end
%
% Construct M^*
%    Diagonal Block
                        jPanel = iPanel;
%%        regular quadrature
                        j = (jPanel - 1)*npt + (1: npt);
                        M0((iBody-1)*npts + itarg, (jBody-1)*npts + j) ...
                            = (Kernel(alpha, z(j, jBody), ...
                                      z(itarg, iBody), Nz(j, jBody))...
                                      .*ds(j, jBody).*w)';
                        M0((iBody-1)*npts + itarg, ...
                           (iBody-1)*npts + itarg) = G0Diag(index);
            
%       correction            
                        wCorrect = wCorrectDiag(index, :)';
                        wCorrect(index) = wCorrectDiag(index, index) ...
                                  + log(0.5*ds(i(index), iBody)*dt);
                        GL = GLog(alpha, z(j, jBody), z(itarg, iBody), ...
                                  Nz(j, jBody));
                        MS((iBody-1)*npts + itarg, ...
                           (jBody-1)*npts + j) ...
                                  = (GL.*ds(j, jBody).*w.*wCorrect)';
                        MS((iBody-1)*npts + itarg, ...
                           (iBody-1)*npts + itarg) = 0;
            
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
                        M0((iBody-1)*npts + itarg, (jBody-1)*npts + j) ...
                            = (Kernel(alpha, z(j, jBody), ...
                                      z(itarg, iBody), Nz(j, jBody))...
                                      .*ds(j, jBody).*w)';
            
%        correction            
                        if abs(t(itarg) - tmid) < ActivateProdFac*dt
                            GL = GLog(alpha, z(j, jBody), ...
                                      z(itarg, iBody), Nz(j, jBody));
                            MS((iBody-1)*npts + itarg, ...
                               (jBody-1)*npts + j) ...
                               = (GL.*ds(j, jBody)...
                                    .*w.*wCorrectSub(index,:)')';
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
                        M0((iBody-1)*npts + itarg, ...
                           (jBody-1)*npts + j) ...
                              = (Kernel(alpha, z(j, jBody), ...
                                        z(itarg, iBody), Nz(j, jBody)) ...
                                        .*ds(j, jBody).*w)';
            
%        correction            
                        if abs(t(itarg) - tmid) < ActivateProdFac*dt
                            GL = GLog(alpha, z(j, jBody), ...
                                      z(itarg, iBody), Nz(j, jBody));
                            MS((iBody-1)*npts + itarg, ...
                               (jBody-1)*npts + j) ...
                             = (GL.*ds(j, jBody).*w...
                                  .*wCorrectSup(index,:)')';
                        end
                    end
               else
                    for index = 1: npt
                        itarg = (iBody-1)*(iPanel - 1)*npt + index;
                        for jPanel = 1: nPanel
                            j = (jPanel - 1)*npt + (1: npt);
                            M0((iBody-1)*npts + itarg, ...
                               (jBody-1)*npts + j) ...
                              = (Kernel(alpha, z(j, jBody), ...
                                        z(itarg, iBody), Nz(j, jBody)) ...
                                        .*ds(j, jBody).*w)';
                        end
                    end
               end
            end
        end
    end
    
    DLP = eye(npts*nBody) + M0 + MS;

end

function [Kern] = Kernel(alpha, zSource, zTarget, Nz)
    dR = zSource - zTarget;
    gradU = besselk(1, abs(dR)/alpha).*dR./abs(alpha*dR);
    Kern = real(gradU.*conj(Nz)) /pi;       
end

function [GL] = GLog(alpha, zSource, zTarget, Nz)
    dR = zSource - zTarget;
    GL = besseli(1, abs(dR)/alpha)...
          .*real(dR.*conj(Nz))./(pi*alpha*abs(dR));
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
  



