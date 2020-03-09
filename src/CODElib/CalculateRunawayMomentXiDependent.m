function moment = CalculateRunawayMomentXiDependent(Nxi,Ny,f,xiIntMat,weightsMat)
%CALCULATERUNAWAYMOMENTXIDEPENDENT calculates the runaway
    fMat = reshape(f,Ny,Nxi)';
    ints = sum( (fMat.*xiIntMat) * weightsMat, 2 );
    plusMin = ones(1,Nxi)';
    plusMin(2:2:end) = -1;
    moment = 0.5*sum(ints.*plusMin);
end
