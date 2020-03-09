    function sigma = CalculateSigmaIntegral(gamma,gamma_m)
        % Calculates the integral
        % int_{gamma_m}^{gamma+1-gamma_m} \Sigma(gamma,gamma_1) dgamma_1
        % analytically. This calculation is needed for the full, conservative
        % close-collision operator. Gamma can be a vector or array, gamma_m
        % must be a scalar. If gamma+1-gamma_m< gamma_m, the integral is
        % set to be zero.

        %Handle invalid elements here, so that we don't have to think about the
        %boundaries of integration when calling the function. This is done
        %so integration is possible. Sigma elements are then set to 0.
        ids = gamma<(2*gamma_m-1);
        gamma(ids) = realmax;

        I2 = 2./(gamma-1) .* log( (gamma-gamma_m)./( gamma_m-1 ) );
        I1 = 2./(gamma-1).^2 .*( 1./(gamma_m-1) - 1./(gamma-gamma_m) + I2  );
        I3 = gamma+1-2*gamma_m;

        sigma = gamma.^2 .*(gamma-1)./(gamma+1).*I1 - ...
            (2*gamma.^2 + 2*gamma - 1)./(gamma.^2-1) .*I2 + I3./(gamma.^2-1);

        sigma(ids) = 0;

        % %%% benchmark with numerical integration
        % if length(gamma)==1
        % dsigma = @(gamma1) gamma^2./((gamma^2-1).*(gamma1-1).^2.*(gamma-gamma1).^2)...
        %     .*( (gamma-1)^2 - (gamma1-1).*(gamma-gamma1)./gamma^2 .*...
        %     ( 2*gamma^2 + 2*gamma - 1 - (gamma1-1).*(gamma-gamma1)));
        % sigma2 = integral(dsigma,a,b)
        % end

    end
