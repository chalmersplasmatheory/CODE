function m = MapToXiIntMat(cumintMat,xiGrid,EOverEc,y2,deltaRef)
    %Assumes that xiGrid is uniformly spaced, starting at 0!

    xics = 1/EOverEc * (1./(deltaRef^2*y2) + 1);
    xics(xics>1) = 1; %Remove points close to the bulk that give unphysical results. They are not in the runaway region, anyway, and will not contribute in the end
    dxi = xiGrid(2)-xiGrid(1);
    xiIds = round(xics/dxi)+1; %Find the id of the nearest point (assumes uniformly spaced xiGrid starting at 0!)
    m = cumintMat(:,xiIds);

    %Each column corresponds to a specific y value. Set those for y<y_c
    %to 0 so that they don't contribute to the integral
    yc = 1/(deltaRef*sqrt(abs(EOverEc)-1));
    id = find(y2>yc^2,1);
    m(:,1:id) = 0;        
end