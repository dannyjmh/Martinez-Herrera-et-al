function [diffExp_Obs] = cossgsea2(TotalGenes, ExpressedGenes, TotalTargets, ExpressedTargets)
prop = ExpressedGenes/TotalGenes;

ExpectedTargets = prop*TotalTargets;
diffExp_Obs = ExpectedTargets - ExpressedTargets;
end

