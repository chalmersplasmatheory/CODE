classdef MomentumGrid < handle
    %Class contataining information about the grid, used so that the
    %grid parameters can be "static" in some sense for all different
    %objects, reducing memory and tidiuosness

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%        Properties                %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % %%%%%%%%%%%%%%%%%%%%%%%%
    %  Reference values
    % %%%%%%%%%%%%%%%%%%%%%%%%%
    % reference - Reference object
    
    properties (SetAccess = protected)
        reference
    end
    % %%%%%%%%%%%%%%%%%%%%%%%%
    %  Resolution parameters:
    % %%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Nxi - number of Legendre modes used to represent the pitch-angle
    %       coordinate (determines resolution in xi=cos(theta))
    % Ny - number of points in the grid used to represent momentum
    % yMax - maximum of the momentum grid, in units of y=gamma*v/v_{th}.
    %        Together, Ny and yMax determine the resolution in momentum
    % y - momentum grid p/p_thermal, where p_thermal is treated
    %       classically
    % y2 - y.^2
    % y4 - y.^4
    % x - v/v_thermal
    % x2 = x.^2
    % gamma - relativistic gamma function
    % xy2 = x.*y2
    % dxdy = derivative of x with respect to y
    % d2xdy2 = second derivative of x to y
    % energyMoment = y2.*(gamma-1);

    properties (SetAccess = protected)
        Nxi = 3;%Number of legandre modes to use
        Ny = 5; %Number of y grid points to use
        yMax = 1 % Maximum normalized momentum
        matrixSize %Size of matrix, namely square Nxi*Ny
        y = 1;%Normalized momentum vector
        y2 %y.^2
        y4 %y.^4
        x %Normalized velocity vector
        x2 %x.^2
        gamma %Relativistic gamma function
        xy2 % x.*y2
        dxdy %Derivative of x with respect to y 
        d2xdy2 % Second derivative of x with respect to y
        energyMoment %
    end


    % %%%%%%%%%%%%%%%%%%%%%%%%
    %  Differential Calculation Params: 
    % %%%%%%%%%%%%%%%%%%%%%%%%
    % ddy - differential operator, differenting with respect to y, requiers
    %       vector of length same as y
    % d2dy2 - same as ddy but second derivetive
    % yWeights - integration operator with respect to y, same limitation as
    %            ddy (row vector)
    % ddx - differential operator with respect to x
    % d2dx2 - second order differential operator with respect to x

    properties (SetAccess = protected)
        %Differantial Calculation Params
        ddy %Differential operator with respect to normalized momentum y
        d2dy2 %Second differential operator with respect to normalized momentum y
        yWeights %Integration operator with respect to y
        ddx %Differential operator with respect to normalized velocity x 
        d2dx2 %Second differential operator with respect to normalized velocity x
    end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Numerical grid parameters:
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % yGridMode - Specifies how the points on the momentum grid are chosen
    %       0 = uniform grid on [0, dy*(Ny-1)]
    %       1 = nonuniform grid, remapped using y = scaleFactor*tan(pi*a*s/2)
    %       2 = nonuniform grid, remapped using y = -scaleFactor*ln(a-s)
    %       3 = nonuniform grid, remapped using y = scaleFactor*s/(a-s)
    %       4 = nonuniform grid, remapped using y = s^2 + gridParameter*s,
    %                   where gridParameter>0 (default)
    %       5 = nonuniform grid, smooth tanh step between dense bulk and
    %                   sparse tail. gridParameter controls the spacing in
    %                   the bulk (and is in units of y). gridStepWidth
    %                   controls the width of the smooth step between bulk
    %                   and tail, and gridStepPosition its position (in
    %                   units of y_c). gridParameter=1/15,
    %                   gridStepWidth=1/50 and gridStepPosition=2 is a good
    %                   starting point, not implemented yet
    %       Above, s is a uniform grid and a is a constant close to 1.
    %       When running a CODE+GO calculation, yGridMode must be
    %       either 0 or 4. The reason is that only these two remappings
    %       are able to handle varying grids in the required way.
    % gridParameter - for use with yGridMode 4 and 5
    % gridStepWidth - for use with yGridMode 5
    % gridStepPosition - for use with yGridMode 5
    % yMaxBoundaryCondition - Boundary condition at yMax
    %       1 = Dirichlet: F=0 (default)
    %       2 = Robin: dF/dy + (2/y)*F=0, which forces F to behave like 1/y^2
    %       3 = Do not apply a boundary condition at yMax
    %       4 = Dirichlet: F=0, with a bit of extra d2dy2 added at the last
    %               few grid points to eliminate ringing
    % artificialDissipationStrength - only used with yMaxBoundaryCondition=4
    %                                 to control the amount of ringing
    % artificialDissipationWidth - in momentum (how close to yMax). Only
    %                              used with yMaxBoundaryCondition=4 to
    %                              control the amount of ringing
    %
    properties (SetAccess = private)
        %Numerical grid parameters
        yGridMode = 4; %Switch determening the grid space to use
        gridParameter = 0.03;%in grid mode 4: s.^2+gridParamter*s, in grid mode 5: not documented yet, inGridMode 6:
        gridStepWidth = 1/50;%in grid mode 5: not documented yet
        gridStepPosition = 2;%in grid mode 5: not documented yet
        yMaxBoundaryCondition = 0; %Switch for what boundary condition to use
        artificialDissipationStrength = 1e-2; %Used to remove ringing in yMaxBoundaryConiditon = 4
        artificialDissipationWidth = 1e-1; %Used to remove ringing in yMaxBoundaryConiditon = 4

    end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Smart grid mode guessing params
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ExpectedEOverEc - Used in gridMode 7 to get a rough estimate of bulk and tail, to know where maximum resolution is needed. (Default 0)

    properties
        ExpectedEOverEc = 0;
    end

    properties (Constant)
       VERSION = 1.0
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%        Constructor               %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %creates new momentum grid
    methods
        function this = MomentumGrid(reference,varargin)
            %Constructor, constructs a new momentum grid with a reference
            %object. In case of more arguments, the relevant properties are
            %also set according to standard matlab syntax: ('name',value)
            for k = 1:2:(nargin-2)
                if isa(this.(varargin{k}),'Reference')
                    error('This function does not yet support an update of the reference object.')
                end
                this.(varargin{k}) = varargin{k+1};
            end
            this.reference = reference;
            reference.setMomentumGrid(this);
            this.initializeMomentumGrid();
        end

    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%   momentum grid  %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%       (set)      %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods    

        function setResolution(this,varargin)
            % Sets resolution parameters, standard matlab
            % syntax ('name',value). 
            for k = 1:2:(nargin-1)
                if isa(this.(varargin{k}),'Reference')
                    error('This function does not yet support an update of the reference object.')
                end
                this.(varargin{k}) = varargin{k+1};
            end
            this.initializeMomentumGrid();
        end
        
        function setNxi(this, Nxi)
            this.Nxi = Nxi;
            this.initializeMomentumGrid();
        end

        function setNy(this, Ny)
            this.Ny = Ny;
            this.initializeMomentumGrid();
        end

        function setyMax(this, yMax)
            this.yMax = yMax;
            this.initializeMomentumGrid();
        end

        function setyGridMode(this, yGridMode)
            this.yGridMode  = yGridMode;
            this.initializeMomentumGrid();
        end

        function setgridParameter(this, gridParameter)
            this.gridParameter  = gridParameter;
            this.initializeMomentumGrid();
        end

        function setgridStepWidth(this, gridStepWidth)
            this.gridStepWidth  = gridStepWidth;
            this.initializeMomentumGrid();
        end

        function setgridStepPosition(this, gridStepPosition)
            this.gridStepPosition  = gridStepPosition;
            this.initializeMomentumGrid();
        end

        function setyMaxBoundaryCondition(this, yMaxBoundaryCondition)
            this.yMaxBoundaryCondition  = yMaxBoundaryCondition;
            this.initializeMomentumGrid();
        end

        function setartificialDissipationStrength(this, artificialDissipationStrength)
            this.artificialDissipationStrength  = artificialDissipationStrength;
            this.initializeMomentumGrid();
        end

        function setartificialDissipationWidth(this, artificialDissipationWidth)
            this.artificialDissipationWidth  = artificialDissipationWidth;
            this.initializeMomentumGrid();
        end
    end
    
    methods
        
        function initializeMomentumGrid(this)
            %Initializes the momentum grid, with the desired
            %non-uniformity. Grid modes 0 and 4 are also able to create a
            %grid with the same spacing as another grid, but extending further
            %in momentum. This can be used to effectively "add" points to a
            %grid without affecting the already existing points, thereby
            %avoiding interpolation of an existing distribution.
            %    o.initializeMomentumGrid();
            normIndex = this.Ny;
            yMin = 0;
            
            switch this.yGridMode  % Generate differentiation matrices.
                % The yWeights vector could be multiplied by any vector of function
                % values on the y grid to integrate in y. This is used when computing
                % the runaway/current densities.
                case 0
                    % Uniform grid, high order finite difference
                    scheme = 12;
                    [this.y, this.yWeights, this.ddy, this.d2dy2] = m20121125_04_DifferentiationMatricesForUniformGrid(this.Ny, yMin, (this.Ny-1)/(normIndex-1)*this.yMax, scheme);
                    %                 [y, yWeights, ddy, d2dy2] = m20120222_02_uniformGridWeightsAnd4thOrderFDDifferentiation(Ny, yMin, (Ny-1)/(normIndex-1)*this.yMax);
                case 1
                    % Uniform grid in a remapped variable s, with y = scaleFactor * tan(pi* s/2),
                    % where scaleFactor is automatically chosen so the spacing of
                    % grid points near y=0 is dy.
                    scheme = 12;
                    [s, sWeights, dds, d2ds2] = m20121125_04_DifferentiationMatricesForUniformGrid(this.Ny, 0, 1, scheme);


                    b = 0.98; %Introduce a parameter close to (but not quite) 1, to avoid Inf's and NaN's
                    this.y = tan(pi* b*s/2);
                    dyds = pi*b./(2*(cos(pi*b*s/2)).^2);
                    d2yds2 = 2*((pi*b/2)^2)*(cos(pi*b*s/2).^(-3)).*sin(pi*b*s/2);
                case 2
                    scheme = 12;
                    [s, sWeights, dds, d2ds2] = m20121125_04_DifferentiationMatricesForUniformGrid(this.Ny, 0, 1, scheme);


                    a = 1.006; %Introduce a parameter close to (but not quite) 1, to avoid Inf's and NaN's            
                    this.y = -log(a-s);
                    dyds = 1./(a-s);
                    d2yds2 = 1./(a-s).^2;
                case 3
                    % Uniform grid in a remapped variable s, with y = scaleFactor*s/(a-s),
                    % where scaleFactor is automatically chosen so the spacing of
                    % grid points near y=0 is dy.
                    scheme = 12;
                    [s, sWeights, dds, d2ds2] = m20121125_04_DifferentiationMatricesForUniformGrid(this.Ny, 0, 1, scheme);


                    a=1.13; %Introduce a parameter close to (but not quite) 1, to avoid Inf's and NaN's
                    this.y = s./(a-s);
                    dyds = a./(a-s).^2;
                    d2yds2 = 2*a./(a-s).^3;
                case 4
                    % Uniform grid in a remapped variable s, with y = s^2+gridParameter*s
                    % This remapping can handle flexible grids in the sense
                    % that new grid points can be added, without the need to
                    % interpolate the distribution.
                    scheme = 12;
                    [s, sWeights, dds, d2ds2] = m20121125_04_DifferentiationMatricesForUniformGrid(this.Ny, 0, (this.Ny-1)/(normIndex-1), scheme);

                    if this.gridParameter<=0
                        error('gridParameter must be positive.')
                    end
                    this.y = s.*s + this.gridParameter*s;
                    dyds = 2*s + this.gridParameter;
                    d2yds2 = 2*ones(size(s));

                 case 6
                    % piece-wise defined momentum grid,
                    % for y>y0, dy \propto y^r,
                    % while for y<y0 it takes a cubic shape with continuous
                    % y, dy and d2y across the boundary y=y0.

                    y0 = this.gridStepPosition/this.reference.deltaRef; 
                    if this.gridStepPosition == 0
                       %if gridStepPosition = 0, set to near critical momentum
                       if max(this.ExpectedEOverEc) > 1
                           pc = 1/sqrt(this.ExpectedEOverEc-1);
                           y0 = 0.9*this.reference.deltaRef^(-1)*pc;
                       else
                           y0 = 5;
                       end
                    end

                    if y0>this.yMax
                       y0 = 0.9*this.yMax;
                    end

                    scheme = 12;                
                    [s, sWeights, dds, d2ds2] = m20121125_04_DifferentiationMatricesForUniformGrid(this.Ny, 0, (this.Ny-1)/(normIndex-1), scheme);

                    r  = this.gridParameter; % a good value that turns out to be fairly robust is r=1.5 
                                       % Note that I hijacked the meaning of gridParameter which is here a large number 
                    s0 = this.gridStepWidth; % fraction of Ny spent on y<y0

                    if r == 1
                       error('gridParameter = 1 not supported (the momentum grid mapping becomes singular). Choose larger or smaller.')
                    end

                    % Calculate the value s0_min which yields the highest grid
                    % density at y=0 and s0_crit above which y is non-monotonic
                    eta = (1-(y0/this.yMax)^(r-1))/(r-1);
                    s0_min = fzero(@(s) (r*eta^2 + 4*eta+6)*s^3+(r*eta^2-4*eta-18)*s^2+18*s-6,0.5);
                    if r>= 2/3
                       k_crit = 2/r*(1+(3*r/2-1)^(1/3));
                    else
                       k_crit = 2/r*(1-sqrt(1-3*r/2));
                    end
                    s0_crit = k_crit/(eta+k_crit);
                    if (s0_min > 0) && (s0_min < s0_crit) 
                       s0_new = s0_min;
                    else
                       s0_new = 0.99*s0_crit;
                    end

                    if s0 > s0_new
                       % if s0 exceeds s0_min or s0_crit, 
                       % set equal to the lower one
                       warning('gridStepWidth is set too large (yielding non-monotonic momentum grid). Setting to maximum value...')
                       s0 = s0_new;
                    end

                    k = eta*s0/(1-s0);

                    rho = 0.01;
                    kappa = (y0/this.yMax)^r;

                    aaa = (6+r*k^2-4*k)/2;
                    bbb = 3*k-3-r*k^2;
                    ccc = (r*k^2+2-2*k)/2;

                    ySBelowS0      = y0     *(aaa*s/s0 +   bbb*(s/s0).^2 +   ccc*(s/s0).^3);
                    dydsSBelowS0   = y0/s0  *(aaa      + 2*bbb*s/s0      + 3*ccc*(s/s0).^2);
                    d2yds2SBelowS0 = y0/s0^2*(           2*bbb           + 6*ccc*s/s0     );

                    ySAboveS0      = y0./( 1 - k*(r-1)*(s/s0-1) ).^(1/(r-1));   % = y0./( 1-(1-(y0/yMax)^(r-1))*(s-s0)./(1-s0)).^(1/(r-1));
                    dydsSAboveS0   = y0/s0*k * (ySAboveS0/y0).^r; % = y0/((r-1)*(1-s0))*(ySAboveS0/y0).^r*(1-(y0/yMax)^(r-1));
                    d2yds2SAboveS0 = r./ySAboveS0 .* dydsSAboveS0.^2;

                    this.y      = ySBelowS0.* (s<=s0) + ySAboveS0 .*(s>s0);
                    dyds   = dydsSBelowS0.*(s<=s0) +dydsSAboveS0.*(s>s0);
                    d2yds2 = d2yds2SBelowS0.*(s<=s0) + d2yds2SAboveS0.*(s>s0);
                    this.y(end) = this.yMax;
                otherwise
                    error('Invalid yGridMode.')
            end
            if this.yGridMode > 0
                %Additional remapping of variables (common for non-uniform grids)

                this.ddy = diag(1./dyds)*dds;
                this.d2dy2 = -diag(d2yds2 ./ (dyds.^3)) * dds + diag((1./dyds).^2)*d2ds2;
                this.yWeights = dyds .* sWeights;
                this.d2dy2(end,:) = 0;  % These values may be NaN otherwise

                % Next, scale the grid so that the final grid point is at yMax:
                %scaleFactor = dy / (y(2) - y(1));
                %scaleFactor = this.yMax / max(y);

                %When adding grid points to an existing grid, we normalize to
                %the point representing the yMax of the old grid, so that all
                %the grid points remain constant and we avoid the need to
                %interpolate the distribution. Normally, normIndex=Ny, and we
                %get the normal behavior.
                scaleFactor = this.yMax/this.y(normIndex);
                this.y = this.y * scaleFactor;
                this.ddy = this.ddy / scaleFactor;
                this.d2dy2 = this.d2dy2 / (scaleFactor*scaleFactor);
                this.yWeights = this.yWeights * scaleFactor;
                this.yWeights(end)=0;
            end
            % Make sure y and yWeights are row vectors:
            this.y = this.y(:)';
            this.yWeights = this.yWeights(:)';

            %calculate dependent variables.
            this.y2 = this.y.*this.y;
            this.y4 = this.y2.*this.y2;
            this.x = this.y ./ sqrt(1+this.reference.deltaRef ^ 2 *this.y2);
            this.gamma = sqrt(1+this.reference.deltaRef^2*this.y2);
            this.xy2 = this.y2.*this.x;
            this.dxdy = (1+this.reference.deltaRef^2*this.y2).^(-3/2);
            this.d2xdy2 = -3*this.reference.deltaRef^2*this.y .* (1+this.reference.deltaRef^2*this.y2).^(-5/2);
            this.x2=this.x.*this.x;
            this.ddx = diag(1./this.dxdy)*this.ddy;
            this.d2dx2 = -diag(this.d2xdy2 ./ (this.dxdy.^3)) * this.ddy + diag((1./this.dxdy).^2)*this.d2dy2;
            this.energyMoment = this.y2.*(this.gamma-1);
            this.matrixSize = this.Nxi*this.Ny;

            if min(this.y) < 0
                error('Negative momenta. Look over the yGridMode!')
            end
            if any(diff(this.y<0))
                error('Non monotone momentum vector. Look over yGridMode!')
            end
        end
    
    end
    
end
