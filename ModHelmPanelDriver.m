%% ModHelmPanelDriver
% Solves the Dirichlet BVP for the modified Helmholtz equation:
%
%     u - alpha^2 \Delta u = 0 in D
%     u = f on \Gamma,
%
% where D is a bounded, multiply-connected domain with smooth boundary.
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

    alpha = 1.;
    uExact = @(r, alpha) besselk(0, r/alpha);
    zSing = 1.5 + 1.5*1i;   % Location of Singularity; 
                            % Used to generate exact solution
                            % Should be located outside domain.
   
   
%% Construct Domain
%
    nPanel = 80;    % Number of Panels
    npt = 16;       % Number of nodes; should be 16 or 32
    
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

    curveParam(1) = struct('shape', 'star', 'parameters', [1 .2 5], ...
                           'tiltAngle', 0, 'centre', 0);
%     curveParam(1) = struct('shape', 'circle', 'parameters', [1 ], ...
%                            'tiltAngle', 0, 'centre', 0);
    curveParam(2) = struct('shape', 'circle', 'parameters', [.3], ...
                        'tiltAngle', 0, 'centre', 0);
                    
    nBody = length(curveParam);
    nPoints = npt * nPanel * nBody;
    
    isUnbounded = false;
    [t, z, dz, d2z, zP, ds, Nz, kappa] ...
        = buildBoundariesPanel(nPanel, npt, T, curveParam, isUnbounded);
    
%% Display Info 
    disp(['Number of Panels: ', num2str(nPanel)])
    disp(['Number of Bodies: ', num2str(nBody)])
    disp(['Number of Nodes: ', num2str(npt)])
    disp(['Total Unknowns: ', num2str(nPoints)])
    disp(['Domain Shape: ', curveParam(:).shape])
               
% 
% Plot domain
    figure()
    subplot(1, 2, 1)
        hold on
        for kBody = 1: nBody
            plot(real(z(:, kBody)), imag(z(:, kBody)), 'k', 'LineWidth', 2)
            plot(real(zP(:, kBody)), imag(zP(:, kBody)), 'r*')
        end
        xlabel('x')
        ylabel('y')
        title('Domain')
    subplot(1, 2, 2)
        hold on
        for kBody = 1: nBody
            plot(kappa(:, kBody), 'k', 'LineWidth', 2)
        end
        xlabel('i')
        ylabel('\kappa')
        title('Curvature') 
        
%% Build grid
% Embed in Box and add a uniform grid
%
    M = 100;    % Build MxM grid
    [xBox, yBox, igrid, LGammaP] ...
              = buildBoxPanel(M, nPanel, npt, nBody, w, z, dz, ds, kappa);
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
    
    
    [DLP, M0, MS] = dlpYukawaPanelMatrix(alpha, nPanel, npt, nBody,...
                                       t, T, w, W, z, ds, Nz, kappa);
                             
    f = uExact(abs(z - zSing), alpha);
    f = -2*alpha^2*f;
    rho = DLP\reshape(f, [], 1);
    rho = reshape(rho, [], nBody);
    
    figure()
    hold on
    for kBody = 1: nBody
        plot(t, rho(:, kBody))
    end
    xlabel('t')
    ylabel('\rho')
    title('Layer Potential Density')
    
    
%% Error Check 
% Check solution at distant target points; ztarg should be inside domain
% well away from boundary

    ztarg = 0.6 + 0.05*exp(1i*(0:15)*2*pi/16);
    figure(1)
    subplot(1, 2, 1)
        plot(real(ztarg), imag(ztarg), 'bo')
    
    uE = uExact(abs(ztarg - zSing), alpha);
    uMax = max(abs(uE));
    uCalc = zeros(size(ztarg));
    for i = 1: length(ztarg)
        uCalc(i) = dlpYukawaPanelEval(alpha, nPanel, npt, nBody, w, z, ...
                                      Nz, ds, rho, ztarg(i));
    end
    disp(['Relative error on sample points = ', ...
           num2str(max(abs(uCalc - uE)/uMax))])
    disp(' ')
    
%% Evaluate Solution on Grid 
%
% Evaluate at regular grid points
    inDomain = find(igrid==1);
    uBox = zeros(size(zBox));
    for index = inDomain'
        zTarg = zBox(index);
        uBox(index) = dlpYukawaPanelEval(alpha, nPanel, npt, nBody, w, ...
                                         z, Nz, ds, rho, zTarg);
    end

%
% Evaluate at near singular grid points
    
    NearSingular = find(igrid > 1);
    nPtInDomain = numel(inDomain) + numel(NearSingular);
    disp(['Number of grid points in domain = ', ...
          num2str(nPtInDomain)])
    zDisMin = 1.d10;
    for index = NearSingular'
        zTarg = zBox(index);
        zmin = min(abs(z - zTarg));
        zDisMin = min(zDisMin, zmin);
        uBox(index) = dlpYukawaPanelNSEval(alpha, nPanel, npt, nBody, ...
                            w, z, zP, LGammaP, Nz, dz, ds, rho, zTarg, ...
                            igrid(index)-1);
    end
    disp(['Minimum Distance between grid point and boundary: ', ...
          num2str(zDisMin)])
    
    uExactBox = zeros(size(zBox));
    uExactBox(igrid>0) = uExact(abs(zBox(igrid>0) - zSing), alpha);
    errBox = abs(uBox - uExactBox);
    
    disp(['Error on Grid = ', num2str(max(max(errBox)))])
    disp(' ')
    
    figure()
    subplot(1, 2, 1)
    contourf(xBox, yBox, uBox, 20)
    hold on
    plot(real(z), imag(z), 'k', 'LineWidth', 2)
    colorbar
    xlabel('x')
    ylabel('y')
    title('Solution')
    
    subplot(1, 2, 2)
    contourf(xBox, yBox, log10(errBox), 20)
    hold on
    plot(real(z), imag(z), 'k', 'LineWidth', 2)
    colorbar
    caxis([-16 0])
    xlabel('x')
    ylabel('y')
    title('Error')
    
