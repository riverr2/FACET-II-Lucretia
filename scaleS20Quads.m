function PS = scaleS20Quads(PS, quadPS, Escale) 
    for PSele = quadPS
        PS(PSele).Ampl = PS(PSele).Ampl*Escale/10;
    end
end