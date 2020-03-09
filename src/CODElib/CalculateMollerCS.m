function cs = CalculateMollerCS(g,g1)
    %Calculates the Møller cross-section (cs) for large-angle
    %electron-electron collisions between an incoming electron with initial
    %gamma=g1 (relatistic gamma function) and an electron at rest, which
    %acquires gamma = g in the collision. g and g1 can be arrays (of the
    %same size).

    g12 = g1.^2;
    preFactor = g12./((g12-1).*(g-1).^2 .* (g1-g).^2);
    sqBracket = 2*g12 + 2*g1 - 1 - (g-1).*(g1-g);
    cs = preFactor .* ( (g1-1).^2 - (g-1).*(g1-g)./g12 .*sqBracket );
    cs(isnan(cs)) = 0;
end
