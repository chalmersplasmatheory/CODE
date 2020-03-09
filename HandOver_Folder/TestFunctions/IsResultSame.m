function isSame = IsResultSame(Dist1,y1, Dist2,y2,Maxwellian, Maxwelly)
    %ISRESULTSAME Returns true if two distribution functions are the same

    % if the instantenous maxwellian distribution (Maxwellian) given at the
    % points Maxwelly is supplied, then the Disttributions are subracted with
    % the maxwellian and the comparision is done with absolute value, otherwise
    % the relative error is used

    tol = 1e-4;
    relTol = 1e-2;
    if nargin == 6
        Maxwellian1 = interp1(Maxwelly,Maxwellian,y1, 'spline');
        Maxwellian2 = interp1(Maxwelly,Maxwellian,y2, 'spline');

        Dist1 = Dist1-Maxwellian1;
        Dist2 = Dist2-Maxwellian2;

        Dist2 = interp1(y2,Dist2,y1, 'spline');

        error = abs(Dist2-Dist1);
        isSame = max(error)<tol;
    else
        if any(y2~=y1)
            if max(y2)>max(y1)
                Dist2 = interp1(y2,Dist2,y1, 'spline');
            else
                Dist1 = interp1(y1,Dist1,y2, 'spline');
            end
        end
        %dont want to devide by 0, so we define 0/0 as 1 since the error is
        %0, so them divided by each other is 1
        Dist2(Dist1 == 0 & Dist2 == 0) = 1;
        Dist1(Dist1 == 0) = 1;
        Relerror =(Dist2-Dist1)./Dist1;
        error = abs(1-Relerror);
        isSame = max(error)<relTol;
    end
    %%%%%%%%%%%%% for testing  %%%%%%%%%%%%%%%%%%%%
    isSame = max(error);
end

