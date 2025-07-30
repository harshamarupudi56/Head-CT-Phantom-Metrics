function [rdat, xcoord, ycoord, zcoord] = rscan_3D(M0, varargin)
% RSCAN_3D performs a radial scan of a 3D matrix.
% This function averages values along spherical shells centered around the origin.
% The procedure is similar to the 2D rscan but in 3D, where:
% [1] Get coordinates of a sphere around an origin.
% [2] Average values of points where the sphere passes through.
% [3] Change the radius of the sphere and repeat until the radial profile is obtained.

if nargin < 1
    error('Input 3D matrix M0 is required.');
end

% Default parameters
xavg = size(M0,2)/2;
yavg = size(M0,1)/2;
zavg = size(M0,3)/2;  % For 3D, we also consider zavg
Rlim = floor(min(size(M0))/2)-1;
dispFlag = 1;
dispFlagC = 0;
rstep = 1;

% Parse varargin to handle custom input parameters
if exist('varargin','var')
    L = length(varargin);
    if rem(L,2) ~= 0, error('Parameters/Values must come in pairs.'); end
    for ni = 1:2:L
        switch lower(varargin{ni})
            case 'xavg', xavg = varargin{ni+1};
            case 'yavg', yavg = varargin{ni+1};
            case 'zavg', zavg = varargin{ni+1};
            case 'rlim', Rlim = varargin{ni+1};
            case 'dispflag', dispFlag = varargin{ni+1};
            case 'dispflagc', dispFlagC = varargin{ni+1};
            case 'rstep', rstep = varargin{ni+1};
        end
    end
end

% Initialize variables
yxz = size(M0);
Rbnd = floor(min(yxz)/2)-1;
if Rlim > Rbnd, Rlim = Rbnd; end
rdat = zeros(1, floor(Rlim/rstep));  % Pre-allocate the radial data

% Iterate over radial distances
for nRho = 1:rstep:floor(Rlim)
    NOP = round(4*pi*nRho^2);  % Number of points on a spherical shell
    [THETA, PHI] = meshgrid(linspace(0, pi, round(sqrt(NOP))), linspace(0, 2*pi, round(sqrt(NOP))));
    RHO = ones(size(THETA)) * nRho;
    
    % Convert spherical to Cartesian coordinates
    [X, Y, Z] = sph2cart(PHI, pi/2-THETA, RHO);
    X = X + xavg;
    Y = Y + yavg;
    Z = Z + zavg;
    
    % Round to nearest voxel coordinates
    X = round(X);
    Y = round(Y);
    Z = round(Z);
    
    % Ensure coordinates are within bounds of the matrix
    X = max(1, min(X, yxz(2)));
    Y = max(1, min(Y, yxz(1)));
    Z = max(1, min(Z, yxz(3)));
    
    % Calculate the average value for the spherical shell
    dat = M0(sub2ind(yxz, Y, X, Z));  % Get the values from the 3D matrix
    rdat(nRho) = sum(dat) / length(dat);  % Average the values
    
    % Optionally display the progress
    if dispFlag
        h1 = figure(100); clf; box on;
        set(h1,'position',[10 500 400 300]);
        set(h1,'units','pixels');
        set(gca,'units','pixels');
        imagesc(squeeze(M0(:,:,round(zavg)))); axis image; hold on;
        plot3(X, Y, Z, 'y.');
        drawnow;
    end
end

% Output coordinates (useful for plotting or further analysis)
xcoord = X;
ycoord = Y;
zcoord = Z;

% Optionally display the radial profile
if dispFlag
    h2 = figure(101); clf; hold on;
    plot(rdat, 'rx', 'linewidth', 2); axis tight;
    title('Radial Profile');
end
end

