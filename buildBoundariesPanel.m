function [t, z, dz, d2z, zP, ds, Nz, kappa] ...
     = buildBoundariesPanel(nPanel, npt, T, curves, isUnbounded)
% BUILDBOUNDARIESPanel(nPoints, curves, isUnbounded) 
%  Constructs parametrized curves; it is assumed the parametrization is
%  on [0, 2 pi].
%
% NOTE: for simply-connected domains only for now!
% 
% INPUTS:
%   nPanel:
%       number of panels
%   npt:
%       number of nodes per panel 
%   T:
%       canonical Gauss-Legendre nodes
%   curves:
%       a structure array containing information about the parametrization.
%       The number of elements equals the number of component curves  
%       (nBody) inthe boundary. The structure contains the following 
%       fields:
%
%       curves.shape:
%           a string which is 'circle', 'star', 'ellipse', 'star2', 
%           'polygon'
%       curves.parameters:
%           an array containing the following information:
%           for a circle: radius
%           for an ellipse: [a b]
%           for a star: [r dr nArms]
%           for a star2: [ [k]; [ak] ]
%               where r(t) = \sum_k ak cos(k t)
%           for a polygon: [z_1 ... z_n] 
%               where z_1 ~= z_n, listed in counter clockwise direction
%       curve.tiltangle:
%           a rotation angle
%       curve.centre:
%           a complex number representing the centre
%   isUnbounded:
%       true if domain is unbounded, false if bounded (this is a place
%       holder for now)
%
% OUTPUTS:
%   t:
%       parameter values at Gauss-Legendre nodes
%   z, dz, d2z: 
%       a nPoints x nBody matrix consisting of the points along the curve,
%       the first and second derivatives with respect to the 
%       parametrization, respectively. Boundary points are equispaced 
%       in the parameter, and listed in a counter clockwise
%       direction. The orientation of dz is determined by the relation
%       of the curve with respect to the domain.
%   zP:
%       Endpoints of panels on curve
%   ds:
%       |dz/dt|
%   Nz:
%       Outward pointing normal
%   kappa:
%       Curvature
%
%    [~, nBody] = size(curves);  

    nPoints = nPanel*npt;
    dt = 2*pi/nPanel;
    
    t = zeros(npt, nPanel);
    
    z = zeros(npt, nPanel);
    dz = zeros(npt, nPanel);
    d2z = zeros(npt, nPanel);
    zP = zeros(nPanel+1, 1);

    for iPanel = 1: nPanel
        ta = (iPanel-1)*dt;
        t(:, iPanel) = ta + 0.5*dt*(T + 1);
        [z(:, iPanel), dz(:, iPanel), d2z(:, iPanel)] ...
                                    = buildCurve(t(:, iPanel), curves);
        zP(iPanel) = buildCurve(ta, curves);
    end
    zP(nPanel+1) = zP(1);
    t = reshape(t, nPoints, []);
    z = reshape(z, nPoints, []);
    dz = reshape(dz, nPoints, []);
    d2z = reshape(d2z, nPoints, []);
    
    ds = abs(dz);
    Nz = -1i*dz./abs(dz);                   % outward normal
    kappa = -imag(dz.*conj(d2z))./ds.^3;    % curvature

%
% Adjust orientation for interior component curves, and curves in
% unbounded domains
%     if isUnbounded
%         dz = -dz; 
%     elseif nBody > 1
%         for kBody = 2: nBody
%             dz(:, kBody) = -dz(:, kBody);
%         end
%     end
     
end

function [z, dz, d2z] = buildCurve(t, curve)

    z0 = curve.centre;
    rotate = exp(1i*curve.tiltAngle);
    
    switch lower(curve.shape)
        case 'circle'
            r = curve.parameters;
            z = z0 + r*exp(1i*t);
            dz = r*1i*exp(1i*t);
            d2z = -r*exp(1i*t); 
      case 'ellipse'
            a = curve.parameters(1);
            b = curve.parameters(2);
            z = z0 + (a*cos(t) + 1i*b*sin(t))*rotate;
            dz = (-a*sin(t) + 1i*b*cos(t))*rotate;
            d2z = (-a*cos(t) - 1i*b*sin(t))*rotate;
        case 'star'
            a = curve.parameters(1);
            b = curve.parameters(2);
            nArms = curve.parameters(3);
            r = sqrt(a^2 + b^2 + 2*a*b*cos(nArms*t));
            dr = -a*b*nArms*sin(nArms*t)./r;
            d2r = -a*b*(nArms^2)*cos(nArms*t) - dr.^2;
            d2r = d2r./r;
            z = z0 + r.*exp(1i*t)*rotate;
            dz = (dr.*exp(1i*t) + 1i*r.*exp(1i*t))*rotate;
            d2z = (d2r.*exp(1i*t) + 2*1i*dr.*exp(1i*t) ...
                   - r.*exp(1i*t))*rotate;
       case 'star2'
            coeff = curve.parameters;
            k = coeff(1,:);
            a_k = coeff(2, :);
            m = length(a_k);
            r = zeros(size(t));
            dr = zeros(size(t));
            d2r = zeros(size(t));
            for i = 1: m
                r = r + a_k(i)*cos(k(i)*t);
                dr = dr - k(i)*a_k(i)*sin(k(i)*t);
                d2r = d2r -k(i)^2*a_k(i)*cos(k(i)*t);
            end
            z = z0 + r.*exp(1i*t)*rotate;
            dz = (dr.*exp(1i*t) + 1i*r.*exp(1i*t))*rotate;
            d2z = (d2r.*exp(1i*t) + 2*1i*dr.*exp(1i*t) ...
                   - r.*exp(1i*t))*rotate;
       otherwise
            error('Unidentified Curve Type')
    end
end