function res = BC_River_System(Ya, Yb, Config)
    % set boundary conditions for bvp4c
    % Ya = conditions at top (x=0)
    % Yb = conditions at bottom (x=Lbottom)
    
    % odd index: concentration (Y(1) is O2 top, Y(2) is dO2dx top)
    res = [
        % Top boundary (Dirichlet: C = C_init)
        Ya(1)  - Config.O2init;
        Ya(3)  - Config.SO4init;
        Ya(5)  - Config.CH4init;
        Ya(7)  - Config.DICinit;
        Ya(9)  - Config.HCO3init;
        Ya(11) - Config.Calcium;
        Ya(13) - Config.Feinit;
        Ya(15) - Config.HSinit;
        
        % Bottom boundary (Neumann: dC/dx = 0)
        % gradients are the even indices
        Yb(2);
        Yb(4);
        Yb(6);
        Yb(8);
        Yb(10);
        Yb(12);
        Yb(14);
        Yb(16);
    ];
end