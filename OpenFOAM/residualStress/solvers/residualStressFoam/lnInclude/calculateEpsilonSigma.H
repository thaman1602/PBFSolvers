epsilon = symm(gradD);

sigma = 2*mu*epsilon + lambda*(I*tr(epsilon)) - I*threeKalpha*(solidResT-T0);
