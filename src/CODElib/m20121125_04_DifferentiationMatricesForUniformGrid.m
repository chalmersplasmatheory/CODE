function [x, w, D, DD] = m20121125_04_DifferentiationMatricesForUniformGrid(N, xMin, xMax, scheme)
% Finite difference differentiation matrices and integration weights for a
% uniform grid.
%
% Created by Matt Landreman,
% Massachusetts Institute of Technology, Plasma Science & Fusion Center, 2012.
%
% Inputs:
%   N = number of grid points.
%   xMin = minimum value in the domain.
%   xMax = maximum value in the domain.
%   scheme = switch for controlling order of accuracy for differentiation
%            and handling of endpoints.
%
% Options for scheme:
% 0 =  The domain [xMin, xMax] is assumed to be periodic. A 3-point stencil
%      is used everywhere. A grid point will be placed at xMin but not
%      xMax.
% 1 =  Same as scheme=0, except a grid point will be placed at xMax but not
%      xMin.
% 2 =  The domain [xMin, xMax] is assumed to be non-periodic. A 3-point
%      stencil is used everywhere.  The first and last row of the
%      differentiation matrices will use one-sided differences, so they
%      will each have a non-tridiagonal element.
% 3 =  The same as scheme=2, except that the first differentiation matrix
%      will use a 2-point 1-sided stencil for the first and last elements
%      so the matrix is strictly tri-diagonal.  The 2nd derivative matrix
%      is the same as for option 2, since it is not possible to compute
%      the 2nd derivative with only a 2-point stencil.
% 10 = The domain [xMin, xMax] is assumed to be periodic. A 5-point stencil
%      is used everywhere. A grid point will be placed at xMin but not
%      xMax.  This option is like scheme=0 but more accurate.
% 11 = Same as scheme=10, except a grid point will be placed at xMax but
%      not xMin.  This option is like scheme=1 but more accurate.
% 12 = The domain [xMin, xMax] is assumed to be non-periodic. A 5-point
%      stencil is used everywhere.  The first two and last two rows of
%      the differentiation matrices will then each have non-pentadiagonal
%      elements.
% 13 = The same as option 12, except that 3-point stencils are used for the
%      first and last rows of the differentiation matrices, and 4-point
%      stencils are used for the 2nd and penultimate rows of the
%      differentiation matrices.  With this option, both differentiation
%      matrices are strictly penta-diagonal.
% 20 = The domain [xMin, xMax] is assumed to be periodic. Spectral
%      differentiation matrices are returned. A grid point will be placed
%      at xMin but not xMax.
% 21 = Same as scheme=20, except a grid point will be placed at xMax but not
%      xMin.
%
% Outputs:
%   x = column vector with the grid points.
%   w = column vector with the weights for integration using the trapezoid rule.
%   D = matrix for differentiation.
%   DD = matrix for the 2nd derivative.

% Validate input:
if ~isfinite(N)
    error('N must be finite')
end
if ~isfinite(xMin)
    error('xMin must be finite')
end
if ~isfinite(xMax)
    error('xMax must be finite')
end
if numel(N) ~= 1
    error('N must be a single integer')
end
if numel(xMin) ~= 1
    error('xMin must be a single number')
end
if numel(xMax) ~= 1
    error('xMax must be a single number')
end
if numel(scheme) ~= 1
    error('scheme must be a single integer')
end
if xMax<xMin
    error('xMax must be larger than xMin.')
end
if xMax==xMin
    error('xMax cannot equal xMin.')
end
if N<2
    error('N must be at least 2.')
end
if N ~= round(N)
    error('N must be an integer.')
end

% Create grid:
switch scheme
    case {0, 1, 10, 11, 20, 21}
        % Periodic
        x=linspace(xMin, xMax, N+1)';
    otherwise
        % Aperiodic
        x=linspace(xMin, xMax, N)';
end
switch scheme
    case {0, 10, 20}
        % Drop point at xMax:
        x(end)=[];
    case {1, 11, 21}
        % Drop point at xMin:
        x(1)=[];
end
dx=x(2)-x(1);
dx2=dx*dx;

% Create integration weights for the trapezoid rule:
w=ones(size(x));
switch scheme
    case {2,3,12,13}
        % Domain is aperiodic
        w(1)=0.5;
        w(end)=0.5;
end
w=w*dx;

% Set the interior points of the differentiation matrices:
switch scheme
    case {0,1,2,3}
        % 2nd order (3-point stencil):
        if N<3
            error('N must be at least 3 for the 3-point stencil methods')
        end
        D=(diag(ones(N-1,1),1)-diag(ones(N-1,1),-1))/(2*dx);
        DD=(-2*diag(ones(N,1),0) + diag(ones(N-1,1),1)+diag(ones(N-1,1),-1))/dx2;

    case {10,11,12,13}
        % 4th order (5-point stencil):
        if N<5
            error('N must be at least 5 for the 5-point stencil methods')
        end
        D=(-diag((1/6)*ones(N-2,1),2) + diag((4/3)*ones(N-1,1),1) - diag((4/3)*ones(N-1,1),-1) + diag((1/6)*ones(N-2,1),-2))/(2*dx);
        DD=(-diag((1/12)*ones(N-2,1),2) + diag((4/3)*ones(N-1,1),1) +diag((-5/2)*ones(N,1),0) + diag((4/3)*ones(N-1,1),-1) - diag((1/12)*ones(N-2,1),-2))/(dx2);

    case {20, 21}
        % Create spectral differentiation matrices.
        % Here I've cut and pasted code from the fourdif.m routine in the
        % DMSuite package by S.C. Reddy and J.A.C. Weideman, available at
        % http://www.mathworks.com/matlabcentral/fileexchange/29
        % or here:
        % http://dip.sun.ac.za/~weideman/research/differ.html

        h=2*pi/N;                                % grid spacing
        kk=(1:N-1)';
        n1=floor((N-1)/2); n2=ceil((N-1)/2);

        % Build first derivative matrix D:
        if rem(N,2)==0                         % compute first column of 1st derivative matrix
            topc=cot((1:n2)'*h/2);
            col1=[0; 0.5*((-1).^kk).*[topc; -flipud(topc(1:n1))]];
        else
            topc=csc((1:n2)'*h/2);
            col1=[0; 0.5*((-1).^kk).*[topc; flipud(topc(1:n1))]];
        end;
        row1=-col1;                            % first row
        D = 2*pi/(xMax-xMin)*toeplitz(col1,row1);

        % Build second derivative matrix DD:
        if rem(N,2)==0                         % compute first column of 2nd derivative matrix
            topc=csc((1:n2)'*h/2).^2;
            col1=[-pi^2/3/h^2-1/6; -0.5*((-1).^kk).*[topc; flipud(topc(1:n1))]];
        else
            topc=csc((1:n2)'*h/2).*cot((1:n2)'*h/2);
            col1=[-pi^2/3/h^2+1/12; -0.5*((-1).^kk).*[topc; -flipud(topc(1:n1))]];
        end;
        row1=col1;                             % first row
        DD=(2*pi/(xMax-xMin))^2*toeplitz(col1,row1);

    otherwise
        error('Invalid input value for ''scheme''')
end


% Handle endpoints:
switch scheme
    case {0,1}
        % 3-point stencil, periodic
        D(1, end) = -1/(2*dx);
        D(end, 1) = 1/(2*dx);

        DD(1, end) = 1/dx2;
        DD(end, 1) = 1/dx2;

    case 2
        % 3-point stencil, aperiodic
        D(1,1)=-1.5/dx;
        D(1,2)=2/dx;
        D(1,3)=-0.5/dx;

        D(end,end)=1.5/dx;
        D(end,end-1)=-2/dx;
        D(end,end-2)=0.5/dx;

        DD(1,1)=1/dx2;
        DD(1,2)=-2/dx2;
        DD(1,3)=1/dx2;

        DD(end,end)=1/dx2;
        DD(end,end-1)=-2/dx2;
        DD(end,end-2)=1/dx2;

    case 3
        % Aperiodic.
        % 2-point stencil for the first and last rows of the first
        % differentiation matrix, so the matrix is strictly tri-diagonal.
        % The 2nd derivative matrix is the same as for scheme=0 (i.e. not
        % strictly tri-diagonal) since it is not possible to approximate
        % the 2nd derivative with a 2-point stencil.

        D(1,1)=-1/dx;
        D(1,2)=1/dx;

        D(end,end)=1/dx;
        D(end,end-1)=-1/dx;

        DD(1,1)=1/dx2;
        DD(1,2)=-2/dx2;
        DD(1,3)=1/dx2;

        DD(end,end)=1/dx2;
        DD(end,end-1)=-2/dx2;
        DD(end,end-2)=1/dx2;

    case {10,11}
        % 5-point stencil, periodic

        D(1, end) = -(4/3)/(2*dx);
        D(1, end-1) = (1/6)/(2*dx);
        D(2, end) = (1/6)/(2*dx);

        D(end, 1) = (4/3)/(2*dx);
        D(end, 2) = -(1/6)/(2*dx);
        D(end-1, 1) = -(1/6)/(2*dx);

        DD(1, end) = (4/3)/dx2;
        DD(1, end-1) = -(1/12)/dx2;
        DD(2, end) = -(1/12)/dx2;

        DD(end, 1) = (4/3)/dx2;
        DD(end, 2) = -(1/12)/dx2;
        DD(end-1, 1) = -(1/12)/dx2;

    case 12
        % 5-point stencil, aperiodic

        D(1,1)= -25/(12*dx);
        D(1,2)= 4/(dx);
        D(1,3)=-3/dx;
        D(1,4)=4/(3*dx);
        D(1,5)=-1/(4*dx);

        D(2,1)= -1/(4*dx);
        D(2,2)= -5/(6*dx);
        D(2,3)=3/(2*dx);
        D(2,4)=-1/(2*dx);
        D(2,5)=1/(12*dx);

        D(end,end)= 25/(12*dx);
        D(end,end-1)= -4/(dx);
        D(end,end-2)=3/dx;
        D(end,end-3)=-4/(3*dx);
        D(end,end-4)=1/(4*dx);

        D(end-1,end)= 1/(4*dx);
        D(end-1,end-1)= 5/(6*dx);
        D(end-1,end-2)=-3/(2*dx);
        D(end-1,end-3)=1/(2*dx);
        D(end-1,end-4)=-1/(12*dx);


        DD(1,1)=35/(12*dx2);
        DD(1,2)=-26/(3*dx2);
        DD(1,3)=19/(2*dx2);
        DD(1,4)=-14/(3*dx2);
        DD(1,5)=11/(12*dx2);

        DD(2,1)=11/(12*dx2);
        DD(2,2)=-5/(3*dx2);
        DD(2,3)=1/(2*dx2);
        DD(2,4)=1/(3*dx2);
        DD(2,5)=-1/(12*dx2);

        DD(end,end)=35/(12*dx2);
        DD(end,end-1)=-26/(3*dx2);
        DD(end,end-2)=19/(2*dx2);
        DD(end,end-3)=-14/(3*dx2);
        DD(end,end-4)=11/(12*dx2);

        DD(end-1,end-0)=11/(12*dx2);
        DD(end-1,end-1)=-5/(3*dx2);
        DD(end-1,end-2)=1/(2*dx2);
        DD(end-1,end-3)=1/(3*dx2);
        DD(end-1,end-4)=-1/(12*dx2);

    case 13
        % Aperiodic.
        % 3-point stencil for the first and last rows of the
        % differentiation matrices, and 4-point stencil for the 2nd and
        % penultimate rows of the differentiation matrices, so the matrices
        % are strictly penta-diagonal.

        D(1,1)=-1.5/dx;
        D(1,2)=2/dx;
        D(1,3)=-0.5/dx;

        D(end,end)=1.5/dx;
        D(end,end-1)=-2/dx;
        D(end,end-2)=0.5/dx;

        D(2,1)=-1/(3*dx);
        D(2,2)=-1/(2*dx);
        D(2,3)=1/(dx);
        D(2,4)=-1/(6*dx);

        D(end-1,end-0)=1/(3*dx);
        D(end-1,end-1)=1/(2*dx);
        D(end-1,end-2)=-1/(dx);
        D(end-1,end-3)=1/(6*dx);

        DD(1,1)=1/dx2;
        DD(1,2)=-2/dx2;
        DD(1,3)=1/dx2;

        DD(end,end)=1/dx2;
        DD(end,end-1)=-2/dx2;
        DD(end,end-2)=1/dx2;

        % It turns out that the 4-point stencil for the second derivative
        % has a weight of 0 for the most distant point, making it identical
        % to the 3-point stencil:

        DD(2,1)=1/(dx2);
        DD(2,2)=-2/(dx2);
        DD(2,3)=1/(dx2);
        DD(2,4)=0;

        DD(end-1,end-0)=1/(dx2);
        DD(end-1,end-1)=-2/(dx2);
        DD(end-1,end-2)=1/(dx2);
        DD(end-1,end-3)=0;

    case {20, 21}
        % Nothing to be done here.

    otherwise
        error('Invalid input value for ''scheme''')
end

end
