function [pValUpReg, pValDnReg, pValDn_Up, ppWT, pmWT, pupWT, pdownWT, UpRegCILow, UpRegCIHigh, DnRegCILow, DnRegCIHigh] = cossgsea(TotalExpected, UpExpected, DownExpected, TotalObserved, UpObserved, DownObserved, textDisp, verbose)

% WT
ppWT= UpExpected/TotalExpected; %cufflink # of genes up / # total of genes
pmWT=DownExpected/TotalExpected; %cufflink # of genes down / # total of genes

% WT
N= TotalObserved; % # total of targets (equivalent to total # of genes in a pathway)
upWT  = UpObserved; % # total of targets that are up 
pupWT =upWT./N;

% WT
downWT=DownObserved; % # total of targets that are down
pdownWT=downWT./N;

if verbose,
  disp([textDisp ' UpReg expected ratio = ' num2str(ppWT) ' DnReg expected ratio = ' num2str(pmWT) ])
  disp([textDisp ' UpReg observed ratio = ' num2str(pupWT) ' DnReg observed ratio = ' num2str(pdownWT) ])
end;

for i=1:1
    % WT
    x=[0:1:N(i)];
    observedupWTPDF=binopdf(x,N(i),pupWT(i));
    expectedupWTPDF=binopdf(x,N(i),(ppWT));
    [diffupWTPDF,diffupWT]=pdfDiff(observedupWTPDF,x,expectedupWTPDF,x);
    % Dice si la diferencia entre lo observado (miRNA) menos lo esperado (mRNA) es
    % significativa
    CDF=cumsum(diffupWTPDF);
    lowerWT=diffupWT(max(find(CDF<0.025)));
    upperWT=diffupWT(max(find(CDF<0.975)));
    if lowerWT>0
        pvalWT=sum(diffupWTPDF(find(diffupWT<=0)));
    elseif upperWT<0
        pvalWT=sum(diffupWTPDF(find(diffupWT>=0)));
    else
        pvalWT=1;
    end
    if verbose,
        disp([textDisp ' UpReg confidence interval 95% =[' num2str(lowerWT) ',' num2str(upperWT) '] pval=' num2str(pvalWT)])
    end;
    UpRegCILow = lowerWT; 
    UpRegCIHigh = upperWT; 
    pValUpReg = pvalWT;
    
    observeddownWTPDF=binopdf(x,N(i),pdownWT(i));
    expecteddownWTPDF=binopdf(x,N(i),(pmWT));
    %[diffdownWTPDF,diffdownWT]=pdfDiff(expecteddownWTPDF,x,observeddownWTPDF,x);
    [diffdownWTPDF,diffdownWT]=pdfDiff(observeddownWTPDF,x,expecteddownWTPDF,x);
    % Dice si la diferencia entre lo observado (miRNA) menos lo esperado (mRNA) es
    % significativa
    CDF=cumsum(diffdownWTPDF);
    lowerWT=diffdownWT(max(find(CDF<0.025)));
    upperWT=diffdownWT(max(find(CDF<0.975)));
    if lowerWT>0
        pvalWT=sum(diffdownWTPDF(find(diffdownWT<=0)));
    elseif upperWT<0
        pvalWT=sum(diffdownWTPDF(find(diffdownWT>=0)));
    else
        pvalWT=1;
    end
    if verbose,
        disp([textDisp ' DownReg confidence interval 95% =[' num2str(lowerWT) ',' num2str(upperWT) '] pval=' num2str(pvalWT)])
    end
    DnRegCILow = lowerWT; 
    DnRegCIHigh = upperWT; 
    pValDnReg = pvalWT;

    %[totalWTPDF,totalWT]=pdfSum(diffupWTPDF,diffupWT,diffdownWTPDF,diffdownWT);
    [totalWTPDF,totalWT]=pdfDiff(diffdownWTPDF,diffdownWT, diffupWTPDF,diffupWT);
    totalWTCDF=cumsum(totalWTPDF);
    lowerWT=totalWT(max(find(totalWTCDF<0.025)));
    upperWT=totalWT(max(find(totalWTCDF<0.975)));
    
    if lowerWT>0
        pvalWT=sum(totalWTPDF(find(totalWT<=0)));
    elseif upperWT<0
        pvalWT=sum(totalWTPDF(find(totalWT>=0)));
    else
        pvalWT=1;
    end
    pValDn_Up = pvalWT;

    if verbose,
%        figure;
        plot(diffupWT,diffupWTPDF,diffdownWT,diffdownWTPDF,totalWT,totalWTPDF)
        legend('UpregulatedDiff: observed-expected','DownregulatedDiff: observed-expected','DownregulatedDiff-UpregulatedDiff')
        disp([textDisp ' DnReg-UpReg Confidence interval 95% =[' num2str(lowerWT) ',' num2str(upperWT) '] pval=' num2str(pvalWT)])
    end;
end

