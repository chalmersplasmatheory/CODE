function outPls = LegendrePolynomials(l,x)
    %LegendrePolynomials calculates P_l(x)
    %
    % Calculates the legendre polynomials P_i(x) for i=0,1,...,l using
    % Bonnet's recursion formula:
    %   (n+1)*P_{n+1}(x) = (2n+1)*x*P_n(x) - n*P_{n-1}(x),
    % where P_0(x) = 1 and P_1(x) = x.    %
    %
    % Usage:
    %       pls = LegendrePolynomials(l,x)
    %
    % l is the (highest) mode number and x is the coordinate, which must be
    % a row vector (not a matrix). pls has the structure
    %   [ P_0(x) ]
    %   [ P_1(x) ]
    %   [   ...  ]
    %   [ P_l(x) ]
    %
    % Written by Adam Stahl 2014

    %Check input
    if ~(isscalar(l) && isnumeric(l) && isreal(l) && l>=0)
        error('Invalid mode number.');
    end
    if size(x,1) > 1
        error('The argument must be presented as a row vector.');
    end
    if any(x<-1) || any(x>1)
        error('The argument is out of the range [-1,1].');
    end

    %Initialize the output and include the first two modes.
    outPls(l+1,numel(x)) = 0;
    outPls(1,:) = ones(size(x));
    if l>=1
        outPls(2,:) = x;
    end

    for n = 1:l-1
        %The index n reflects the n in Bonnet's formula, but the
        %corresponding row in the array is n+1. We start from n=1, since we
        %want to use the formula to compute n+1=2 (which goes on row n+2=3)
        outPls(n+2,:) = (2*n+1)/(n+1) * x.*outPls(n+1,:) - n/(n+1) * outPls(n,:);
    end

end
