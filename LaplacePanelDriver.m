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

   uExact = @(z) real(z.^2);
   
   
%% Construct Domain
%
    nPanel = 40;    % Number of Panels
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

    curveParam(1) = struct('shape', 'circle', 'parameters', [1], ...
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
    M = 100;    % Build MxM grid
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
    
    
    DLP = dlpLaplacePanelMatrix(nPanel, npt, w, z, ds, Nz, kappa);
    f = 2*uExact(z);
    [SLP, M0, MS] = slpLaplacePanelMatrix(nPanel, npt, t, T, w, W, z, ds);
    f = uExact(z);
%
% add in regularization for SLP
    rowAV = 0.5*repmat(w, nPanel, 1).*ds/pi;
    AV = repmat(rowAV', nPoints, 1);
    SLP = SLP + AV;
    
    rho = SLP\f;
    
    figure()
    plot(t, rho)
    xlabel('t')
    ylabel('\rho')
    title('Layer Potential Density')
    
    
%% Error Check 
% Check solution at distant target points; ztarg should be inside domain
% well away from boundary

    ztarg = 0.25*exp(1i*(0:15)*2*pi/16);
    uE = uExact(ztarg);
    uCalc = zeros(size(ztarg));
    for i = 1: length(ztarg)
%         uCalc(i) = dlpLaplacePanelEval(nPanel, npt, w, z, Nz, ds, rho, ...
%                                        ztarg(i));
        uCalc(i) = slpLaplacePanelEval(nPanel, npt, w, z, ds, rho, ...
                                       ztarg(i), AV);
    end
    disp(['Error on sample points = ', num2str(max(abs(uCalc - uE)))])
    disp(' ')
    
%% Evaluate Solution on Grid 
%
% Evaluate at regular grid points
    inDomain = find(igrid==1);
    uBox = zeros(size(zBox));
    for index = inDomain'
        zTarg = zBox(index);
        uBox(index) = slpLaplacePanelEval(nPanel, npt, w, z, ds, ...
                                          rho, zTarg, AV);
    end

%
% Evaluate at near singular grid points
    
    NearSingular = find(igrid==2);
    for index = NearSingular'
        zTarg = zBox(index);
        
        uBox(index) = slpLaplacePanelNSEval(nPanel, npt, w, z, zP, ...
                                     LGammaP, dz, Nz, ds, rho, zTarg, AV);
    end
    
    uExactBox = zeros(size(zBox));
    uExactBox(igrid > 0) = uExact(zBox(igrid > 0));
    errBox = abs(uBox - uExactBox);
    
    disp(['Error on Grid = ', num2str(max(max(errBox)))])
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
    
