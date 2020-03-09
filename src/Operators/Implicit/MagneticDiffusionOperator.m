classdef MagneticDiffusionOperator < ImplicitOperator
    %MAGNETICDIFFUSIONOPERATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %Local transport sink
        deltaBoverB = 0.001;            % Size of radial magnetic field fluctuations
        safetyFactor = 1.0;             % Safety factor (q) needed for diffusion coefficient estimate
        majorRadius = 1.0;              % Major radius (in meters) for Rechester-Rosenbluth diffusion
        radialScaleLength = 0.1;        % Length scale for radial variation of runaway density (in meters)
        ionMassNumber = 1.0;            % Ion mass number (in proton masses)
    end
    
    methods
        function this = MagneticDiffusionOperator(state, eqSettings, varargin)
            this@ImplicitOperator(state,eqSettings,varargin);
        end
        
        function matrixHasChanged = generateOperatorMatrix(this,iteration)
            matrixHasChanged = 1; %(todo smart implementation to not change in every iteration)
            
            
            
            nueeRef     = this.state.reference.nueeRef;
            T           = this.state.physicalParams.T(iteration);
            Ny = this.state.momentumGrid.Ny;
            Nxi = this.state.momentumGrid.Nxi;

            this.predictedNNZ = (Nxi*Nxi) * Ny;
            this.resetSparseCreator();
            
            switch this.eqSettings.collisionOperator
                case {0,2,3,4}
                    if this.eqSettings.useNonRelativisticCollOp
                        takeNonRelLimit = true;
                    else
                        takeNonRelLimit = false;
                    end
                case {1}
                    takeNonRelLimit = true;
                otherwise
                    error('You have chosen an invalid collision operator.');
            end

            if takeNonRelLimit
                y = this.state.momentumGrid.x;
            else
               y = this.state.momentumGrid.y;
            end
            % Calculate coupling coefficients between L and l when including diffusion-like operator
            if this.deltaBoverB > 0.0
                absCoefficients = zeros(Nxi,Nxi);

                absCoefficients(1,1) = 1.0;
                
                % Special case: L=0
                for l=2:2:(Nxi-1)
                    num = l*(2*l+3)*legendreP(l-2,0) - (2*l+1)*legendreP(l,0) - (l+1)*(2*l+1)*legendreP(l+2,0);
                    denom = (2*l-1)*(2*l+1)*(2*l+3);
                    absCoefficients(1,l+1) = 2*num/denom;
                end

                % Special case: l=0
                for L=2:2:(Nxi-1)
                    num = L*(2*L+3)*legendreP(L-2,0) - (2*L+1)*legendreP(L,0) - (L+1)*(2*L+1)*legendreP(L+2,0);
                    denom = (2*L-1)*(2*L+1)*(2*L+3);
                    absCoefficients(L+1,1) = 2*num/denom;
                end

                for L=1:(Nxi-1)
                    for l=1:(Nxi-1)
                        % only nonzero if integrand is even
                        if (mod(l+L,2) == 0)

                            if L == (l+1)
                                res_p1 = 1.0/(2*L+1);
                            else
                                if mod(L,2) == 0
                                    num = (-1)^(0.5*(L+l+2))*factorial(L)*factorial(l+1);
                                    denom = 2^(L+l)*(L-l-1)*(L+l+2)*(factorial(L/2))^2*(factorial(l/2))^2;
                                    res_p1 = num/denom;
                                else
                                    num = (-1)^(0.5*(L+l+2))*factorial(L)*factorial(l+1);
                                    denom = 2^(L+l)*(l-L+1)*(L+l+2)*(factorial((l+1)/2))^2*(factorial((L-1)/2))^2;
                                    res_p1 = num/denom;
                                end
                            end

                            if L == (l-1)
                                res_m1 = 1.0/(2*L+1);
                            else
                                if mod(L,2) == 0
                                    num = (-1)^(0.5*(L+l))*factorial(L)*factorial(l-1);
                                    denom = 2^(L+l-2)*(L-l+1)*(L+l)*(factorial(L/2))^2*(factorial((l-2)/2))^2;
                                    res_m1 = num/denom;
                                else
                                    num = (-1)^(0.5*(L+l))*factorial(L)*factorial(l-1);
                                    denom = 2^(L+l-2)*(l-L-1)*(L+l)*(factorial((l-1)/2))^2*(factorial((L-1)/2))^2;
                                    res_m1 = num/denom;
                                end
                            end

                            absCoefficients(L+1,l+1) = (l+1)*res_p1 + l*res_m1;

                        end
                    end
                end

                for L=0:(Nxi-1)
                    % Including diffusion makes the matrix dense in xi, so we need to make it a special case
                    %%%%% Set parameters for diffusion %%%
                    % Normalized diffusion coefficient
                    meSI = 9.10938356e-31;
                    mpSI = 1.6726219e-27;
                    eSI = 1.60217662e-19;
                    c = 2.997925e8; %m/s
                    
                    veRef = this.state.reference.deltaRef*c; %(todo, added by Albert, check if it is correct units)
                    
                    Dhat = pi*this.safetyFactor*this.majorRadius*(this.deltaBoverB/this.radialScaleLength)^2 * (veRef*0.01/nueeRef);
                    soundSpeed = sqrt( eSI*T / (this.ionMassNumber*mpSI));

                    drOverlB = sqrt(gamma(2:Ny).^2 - 1.0)*this.safetyFactor*c*(meSI/(mpSI*this.ionMassNumber))/(2.5*soundSpeed);

                    Yp = abs(1.0 - 2*drOverlB.^2 + 2.5*drOverlB.^4);
                    Ypp = (1.5*pi*drOverlB).^(-1);

                    Ytot = min(Yp,Ypp);

                    for l=0:(Nxi-1)
                        diffTerm = -Dhat * absCoefficients(L+1,l+1)*Ytot.*(y(2:Ny)./gamma(2:Ny));
                        rowIndices = L*Ny + (2:Ny);
                        colIndices = l*Ny + (2:Ny);

                        this.addToSparse(rowIndices, colIndices, diffTerm);
                    end

                    diffOut(2:Ny) = Dhat * this.radialScaleLength^2 * nueeRef * Ytot.*(y(2:Ny)./gamma(2:Ny));
                    diffOut(1) = 0.0;
                end
            end
            
            this.operatorMatrix = -this.createSparse();
            this.clearSparseResidue();
            
        end
        
    end
    
end

