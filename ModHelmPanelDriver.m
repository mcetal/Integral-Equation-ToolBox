%% ModHelmPanelDriver
% Solves the Dirichlet BVP for the modified Helmholtz equation:
%
%     u - alpha^2 \Delta u = 0 in D
%     u = f on \Gamma,
%
% where D is a bounded, simply-connected domain with smooth boundary.
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
%   dlpYukawaPanelMatrix.m
%   dlpYukawaPanelEval.m
%   dlpYukawaPanelNSEval.m
%
%
% SIMPLY CONNECTED DOMAIN ONLY!
%
   close all; 
   clear;
   clc
   
%% PDE
%

    alpha = 1;
    uExact = @(r, alpha) besselk(0, r/alpha);
    zSing = 1.5 + 1.5*1i;
   
   
%% Construct Domain
%
    nPanel = 10;
    npt = 16;
    nPoints = npt * nPanel;
    dt = 2*pi/nPanel;
    
    if npt == 16
        [T, W] = GaussLegendre16;
    else
        [T, W] = GaussLegendre32;
    end
    
    w = 0.5*dt*W;
    
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

    curveParam(1) = struct('shape', 'circle', 'parameters', [1 ], ...
                        'tiltAngle', 0, 'centre', 0);

    isUnbounded = false;
    [t, z, dz, d2z, zP, ds, Nz, kappa] ...
        = buildBoundariesPanel(nPanel, npt, T, curveParam, isUnbounded);
               
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
    zBox = xBox + 1i*yBox;
    figure()
    mesh(xBox, yBox, igrid)
    hold on
    plot3(real(z), imag(z), 2*ones(size(z)), 'LineWidth', 2)
    
%% Integral Equation
% Build system matrix and solve integral equation
    
    
    [DLP, M0, MS] = dlpYukawaPanelMatrix(alpha, nPanel, npt, t, T, w, W,...
                                         z, ds, Nz, kappa);
                             
                             
    f = uExact(abs(z - zSing), alpha);
    f = -2*alpha^2*f;
    rho = DLP\f;
    
    figure()
    plot(t, rho)
    
    
%%
% Check solution at distant target points; ztarg should be inside domain
% well away from boundary

    ztarg = 0.25*exp(1i*(0:15)*2*pi/16);
    uE = uExact(abs(ztarg - zSing), alpha);
    uCalc = zeros(size(ztarg));
    for i = 1: length(ztarg)
        uCalc(i) = dlpYukawaPanelEval(alpha, nPanel, npt, w, z, Nz, ds, ...
                                      rho, ztarg(i));
    end
    disp(['Error on sample points = ', num2str(max(abs(uCalc - uE)))])
    
%% Evaluate Solution on Grid 
%
% Evaluate at regular grid points
    inDomain = find(igrid==1);
    uBox = zeros(size(zBox));
    for index = inDomain'
        zTarg = zBox(index);
        uBox(index) = dlpYukawaPanelEval(alpha, nPanel, npt, w, z, Nz, ...
                                         ds, rho, zTarg);
    end

%
% Evaluate at near singular grid points
    
    NearSingular = find(igrid==2);
    for index = NearSingular'
        zTarg = zBox(index);
        
        uBox(index) = dlpYukawaPanelNSEval(alpha, nPanel, npt, w,   ...
                            z, zP, LGammaP, Nz, dz, ds, rho, zTarg);
    end
    
    uExactBox = zeros(size(zBox));
    uExactBox(igrid>0) = uExact(abs(zBox(igrid>0) - zSing), alpha);
    errBox = abs(uBox - uExactBox);
    
    disp(['Error on Grid = ', num2str(max(max(errBox)))])
    
    figure()
    contourf(xBox, yBox, log10(errBox), 20)
    colorbar
    caxis([-16 0])
    
