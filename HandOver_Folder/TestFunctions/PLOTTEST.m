for j = 11:12
    load(['TestResults/Test_14-Oct-2019_' num2str(j)])
    for i = 1:length(resCODE.times)
        clf
        figure(1)
        subplot(1,2,1)
        semilogy(resCODE.y,abs(resCODE.f(1:resCODE.Ny,i)),'b')
        hold on
        semilogy(resCODENew.distributions{i}.momentumGrid.y,abs(resCODENew.distributions{i}.f(1:length(resCODENew.distributions{i}.momentumGrid.y))),'--r')
        legend('CODE','CODE Object Oriented')
        subplot(1,2,2)
        plot(1:resCODE.Ny*resCODE.Nxi,...
            (resCODE.f(:,i)...
            -resCODENew.distributions{i}.f(:))...
            ./(resCODE.f(:,i)+(resCODE.f(:,i) == 0)))
        hold on
        plot((1:resCODE.Nxi)*resCODE.Ny,0,'*')
        %find(isinf((resCODE.f(:,i)...
        %    -resCODENew.distributions{i}.f(:))...
        %    ./(resCODE.f(:,i)+(resCODE.f(:,i)==0))))
        norm((resCODE.f(:,i)...
            -resCODENew.distributions{i}.f(:))...
            ./(resCODE.f(:,i)+(resCODE.f(:,i) == 0)))
       pause
    end
end