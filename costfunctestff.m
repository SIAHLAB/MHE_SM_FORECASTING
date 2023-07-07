function [cost] = costfunctestff(mpres,thetf,mwant,mwantd,kw,thet2g2,riwant,rainwd1,rainwd2,rainwd3,rainwd4,rainwd5,rainwd6,tewant,ts,thetantf,alphapost)

cost=norm(mpres-estimodel(thetf,mwant,mwantd,kw,thet2g2,riwant,rainwd1,rainwd2,rainwd3,rainwd4,rainwd5,rainwd6,tewant,ts)')^2 + alphapost*(thetantf-thetf)^2;