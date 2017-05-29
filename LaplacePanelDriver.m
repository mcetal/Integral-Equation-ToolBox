%% LaplacePanelDriver
% Solves the Dirichlet BVP for Laplace's Equation:
%
%     \Delta u = 0 in D
%     u = f on \Gamma,
%
% where D is a bounded, simply-connected domain with smooth boundary.
%
% The integral equation is
%   [I + K] rho = 2 * f
%
% Quadrature is Gauss-Legendre via 
% Helsing & Holst, Variants of an explicit kernel-split panel-based
% Nystrom discretization scheme for Helmholtz boundary value problems,
% (2015)
%
% Represent u using a double layer potential.
%
% Requires the following to run:
%   buildBoundariesPanel.m
%   buildBoxPanel.m
%   dlpLaplacePanelMatrix.m
%   dlpLaplacePanelEval.m
%   dlpLaplacePanelNSEval.m
%   GaussLegendre.m
%
%
% SIMPLY CONNECTED DOMAIN ONLY!
%
   close all; 
   clear;
   clc
   
%% PDE
%

    zstar = 2.5 + 2.4*1i;
    uExact = @(z) real(1./(z - zstar));
    graduExact = @(z) conj(-1./(z-zstar).^2);
%     uExact = @(z) real((z-zstar).^3);
%     graduExact = @(z) conj(3*(z-zstar).^2);
    
% BVP
    iDirichlet = true;
    iSLP = true;   
   
   
%% Construct Domain
%
    nPanel = 100;    % Number of Panels
    npt = 16;       % Number of nodes; should be 16 or 32
    
    nPoints = npt * nPanel;
    dt = 2*pi/nPanel;
    
    [T, W] = GaussLegendre(npt);    
    w = 0.5*dt*W;   % Regular quadrature weights
    
%
% Structure containing information about parametrization
%       curveParam.shape:
%           circle
%           ellipse
%           star
%           star2
%       curveParam.parameters:
%           an array containing the following information:
%           for a circle: radius
%           for an ellipse: [a b]
%           for a star: [r dr nArms]
%           for a star2: [ [k]; [ak] ]

    curveParam(1) = struct('shape', 'circle', 'parameters', [2], ...
                        'tiltAngle', 0, 'centre', 0);

    isUnbounded = false;
    [t, z, dz, d2z, zP, ds, Nz, kappa] ...
        = buildBoundariesPanel(nPanel, npt, T, curveParam, isUnbounded);
    
%% Display Info 
    disp(['Number of Panels: ', num2str(nPanel)])
    disp(['Number of Nodes: ', num2str(npt)])
    disp(['Total Unknowns: ', num2str(nPoints)])
    disp(['Domain Shape: ', curveParam.shape])
               
% 
% Plot domain
    figure()
    subplot(1, 2, 1)
        hold on
        plot(real(z), imag(z), 'k', 'LineWidth', 2)
        hold on
        plot(real(zP), imag(zP), 'r*')
        xlabel('x')
        ylabel('y')
        title('Domain')
    subplot(1, 2, 2)
        hold on
        plot(kappa, 'k', 'LineWidth', 2)
        xlabel('i')
        ylabel('\kappa')
        title('Curvature') 
        
%% Build grid
% Embed in Box and add a uniform grid
%
    M = 200;    % Build MxM grid
    [xBox, yBox, igrid, LGammaP] ...
                 = buildBoxPanel(M, nPanel, npt, w, z, dz, ds);
    disp(['Arc Length of Boundary = ', num2str(sum(LGammaP))])
    disp(' ')

    zBox = xBox + 1i*yBox;
    figure()
    mesh(xBox, yBox, igrid)
    hold on
    plot3(real(z), imag(z), 2*ones(size(z)), 'LineWidth', 2)
    xlabel('x')
    ylabel('y')
    title('Grid Point Flag')
    
%% Integral Equation
% Build system matrix and solve integral equation
    [DH, SH, dsh] = Helsing;
    SH = 0.5*SH;
%
% Set BCs and get System Matrix according to BVP type
    if iDirichlet && ~iSLP
        disp('Dirichlet BVP using DLP')
        f = 2*uExact(z);
     	SysMat = dlpLaplacePanelMatrix(nPanel, npt, w, z, ds, Nz, kappa);
    elseif iDirichlet && iSLP
        disp('Dirichlet BVP using SLP')
        f = uExact(z);
        SysMat = slpLaplacePanelMatrix(nPanel, npt, t, T, w, W, z, ds);
    else
        disp('Neuman BVP')
        f = -2*real(graduExact(z).*conj(Nz));
        SysMat = DslpDnLaplacePanelMatrix(nPanel, npt, w, z, ds, Nz, ...
                                          kappa, false);
    end
    disp(' ')
    
    rho = SysMat \ f;
            
    figure()
    plot(t, rho)
    xlabel('t')
    ylabel('\rho')
    title('Layer Potential Density')
    
    rhoAv = 0;
    for iPanel = 1: nPanel
        i = (iPanel-1)*npt + (1: npt);
        rhoAv = rhoAv + sum(0.5*rho(i).*w.*ds(i)/pi);
    end
    disp(['Average Density (should be zero if Neumann BVP) = ', ...
         num2str(rhoAv)])
    disp(' ')
    
    
%% Error Check 
% Check solution at distant target points; ztarg should be inside domain
% well away from boundary

    ztarg = 0.25*exp(1i*(0:15)*2*pi/16);
    uE = uExact(ztarg);
    uCalc = zeros(size(ztarg));
    for i = 1: length(ztarg)
        if iDirichlet && ~iSLP
            uCalc(i) = dlpLaplacePanelEval(nPanel, npt, w, z, Nz, ds, rho, ...
                                           ztarg(i));
        else
            uCalc(i) = slpLaplacePanelEval(nPanel, npt, w, z, ds, rho, ...
                                           ztarg(i));
        end
    end
    if iDirichlet
        disp(['Error on sample points = ', num2str(max(abs(uCalc - uE)))])
    else
        diff = uCalc - uE;
        err = max(diff) - min(diff);
        disp(['Error on sample points = ', num2str(abs(err))])
    end
    disp(' ')
    
%% Evaluate Solution on Grid 
%
% Evaluate at regular grid points
    inDomain = find(igrid==1);
    uBox = zeros(size(zBox));
    for index = inDomain'
        zTarg = zBox(index);
        if iDirichlet && ~iSLP
            uBox(index) = dlpLaplacePanelEval(nPanel, npt, w, z, Nz, ...
                                              ds, rho, zTarg);
        else
            uBox(index) = slpLaplacePanelEval(nPanel, npt, w, z, ds, ...
                                              rho, zTarg);
        end
    end

%
% Evaluate at near singular grid points
    
    NearSingular = find(igrid==2);
    for index = NearSingular'
        zTarg = zBox(index);
        if iDirichlet && ~iSLP
            uBox(index) = dlpLaplacePanelNSEval(nPanel, npt, w, z, zP, ...
                                      LGammaP, Nz, dz, ds, rho, zTarg);
        else
            uBox(index) = slpLaplacePanelNSEval(nPanel, npt, w, z, zP, ...
                                         LGammaP, dz, Nz, ds, rho, zTarg);
        end
        
    end
    
    uExactBox = zeros(size(zBox));
    uExactBox(igrid > 0) = uExact(zBox(igrid > 0));
    errBox = abs(uBox - uExactBox);
    
    if iDirichlet
        disp(['Error on Grid = ', num2str(max(max(errBox)))])
    else
        err = max(max(errBox(igrid > 0))) - min(min(errBox(igrid > 0)));
        disp(['Error on Grid = ', num2str(abs(err))])
        constant = max(max(errBox));
        errBox(igrid > 0) = errBox(igrid > 0) - constant;
        errBox = abs(errBox);
    end
    disp(' ')
    
    figure()
    contourf(xBox, yBox, log10(errBox), 20)
    hold on
    plot(real(z), imag(z), 'k', 'LineWidth', 2)
    colorbar
    caxis([-16 0])
    xlabel('x')
    ylabel('y')
    title('Error')
    
