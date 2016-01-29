function scriptcossgsea(InFileName, OutFileName, pvalTh, DnTh, UpTh, columns)
    
 %[num,txt,raw] = xlsread(InFileName);
%columns = 837;
%columns = 116;
[num txt] = readmirdata(InFileName,columns);

[len,wid] = size(num);
len = len +1;

IndexesFC = [2];
IndexesmiR=[4:columns];
%IndexesFC = [2 5 8 11 14 17 20 23 26 29];
%IndexesmiR = [32 34 35 36 37];
for j = 1:length(IndexesFC), 
    fileID = fopen(OutFileName,'w');

    % Print header
    indexFC = IndexesFC(j);
    myIndexText = char(txt{indexFC});
    fprintf(fileID,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'miRNA', 'totalExpected', 'expectedUp', 'expectedDn', 'totalObserved', 'observedUp', 'observedDn', 'UpReg expected ratio', 'DnReg expected ratio', 'UpReg observed ratio', 'DnReg observed ratio', 'UpReg confidence interval 95% lower', 'UpReg confidence interval 95% upper', 'Upreg pvalue', 'DnReg confidence interval 95% lower', 'DnReg confidence interval 95% upper', 'Dnreg pvalue', 'Permutation test Up pValue', 'Permutation test Dn pValue');
   
    for k = 1:length(IndexesmiR), 
        
        indexFC = IndexesFC(j); 
        myIndexText = char(txt{indexFC});
    
        indexmiR = IndexesmiR(k); 
        myText = char(txt{indexmiR});

        totalExpected = len-1;
        totalObserved = 0;
        expectedUp = 0;
        expectedDn = 0;
        observedUp = 0;
        observedDn = 0;
    
        %DnTh = -10;
        %UpTh = 10;
        usePval = 1;
        PvalIndex = indexFC+1;  % pvalue
        
        % Uncomment this if you ONLY want to use pvalues
        %if (usePval == 1),
        %    DnTh = 0; UpTh = 0;
        %end;

        for i = 2:len,
            %if (~strcmp(txt(i,indexmiR), ''))
            if (num(i-1,indexmiR-1) == 1)
                totalObserved = totalObserved+1;
                if (num(i-1,indexFC-1) <= DnTh),
                    observedDn = observedDn +1;            
                    if (usePval ==1),
                        if (num(i-1,PvalIndex) > pvalTh),
                            observedDn = observedDn -1;            
                        end;
                    end;
                end;
                if (num(i-1,indexFC-1) > UpTh),
                    observedUp = observedUp +1;            
                    if (usePval ==1),
                        if (num(i-1,PvalIndex) > pvalTh),
                            observedUp = observedUp -1;            
                        end;
                    end;
                end;
            end;
            if (num(i-1,indexFC-1) <= DnTh),
                expectedDn = expectedDn +1;  
                if (usePval ==1),
                    if (num(i-1,PvalIndex) > pvalTh),
                        expectedDn = expectedDn -1;            
                    end;
                end;
            end;
            if (num(i-1,indexFC-1) > UpTh),
                expectedUp = expectedUp +1;            
                if (usePval ==1),
                    if (num(i-1,PvalIndex) > pvalTh),
                        expectedUp = expectedUp -1;            
                    end;
                end;
            end;
        end

        %if (totalObserved == 0),
         %  totalObserved = 1;
         %  observedDn = 1;
        %end;  
 

        disp(['File: ' InFileName]);
        disp(['Fold change method: ' myIndexText]);
        disp(['miRNA: ' myText]);
        fprintf(fileID,'%s',myText);
        disp(['pvalue Threshold: ' num2str(pvalTh)]);
        disp([' ']);
        %disp([myText  ' totalExpected = ' num2str(totalExpected) ' expectedUp = ' num2str(expectedUp) ' expectedDn = ' num2str(expectedDn)]);
        fprintf(fileID,'\t%f\t%f\t%f',totalExpected, expectedUp, expectedDn);
        %disp([myText  ' totalObserved = ' num2str(totalObserved) ' observedUp = ' num2str(observedUp) ' observedDn = ' num2str(observedDn)]);
        fprintf(fileID,'\t%f\t%f\t%f',totalObserved, observedUp, observedDn);
        [pValUpReg, pValDnReg, pValDn_Up, ppWT, pmWT, pupWT, pdownWT, UpRegCILow, UpRegCIHigh, DnRegCILow, DnRegCIHigh] = cossgsea(totalExpected, expectedUp,expectedDn, totalObserved, observedUp, observedDn, myText, 0);
        fprintf(fileID,'\t%f\t%f\t%f\t%f\t%f\t%f\t%e\t%f\t%f\t%e',ppWT, pmWT, pupWT, pdownWT, UpRegCILow, UpRegCIHigh, pValUpReg, DnRegCILow, DnRegCIHigh, pValDnReg);
        %disp(['miR217  '  'pValUpReg = ' num2str(pValUpReg) ' pValDnReg = ' num2str(pValDnReg) ' pValDn-Up = ' num2str(pValDn_Up)]);
        
        % Calculate Fold Enrichment (Immunity 27, 847–859, December 2007)
        %FEup = log(observedUp/totalObserved);
        %FEdn = log(observedDn/totalObserved);
        %disp([' FEup = ' num2str(FEup) ' FEdn = ' num2str(FEdn)]);
        

        % End of Fold Enrichment

        % CossGsea2 (random test)

        [realDifference] = cossgsea2(totalExpected, expectedUp, totalObserved, observedUp);
        [realDifferenceDn] = cossgsea2(totalExpected, expectedDn, totalObserved, observedDn);
        diffExp_ObsUp = zeros(1,10000);
        diffExp_ObsDn = zeros(1,10000);
        for rndI = 1:10000,
            % Now a random vector
            %myRandIdx = unidrnd(totalExpected,1,totalObserved);
            myRandIdx = randperm(totalExpected);
            myRandIdx = myRandIdx(1:totalObserved);
            observedUp = 0;
            observedDn = 0;

           for i = 1:totalObserved,
                if (num(myRandIdx(1,i),indexFC-1) <= DnTh),
                    observedDn = observedDn +1;            
                    if (usePval ==1),
                        if (num(myRandIdx(1,i),PvalIndex) > pvalTh),
                            observedDn = observedDn -1;            
                        end;
                    end;
                end;
                if (num(myRandIdx(1,i),indexFC-1) > UpTh),
                    observedUp = observedUp +1;            
                    if (usePval ==1),
                        if (num(myRandIdx(1,i),PvalIndex) > pvalTh),
                            observedUp = observedUp -1;            
                        end;
                    end;
                end;
           end
               
          [diffExp_ObsUp(rndI)] = cossgsea2(totalExpected, expectedUp, totalObserved, observedUp);
          [diffExp_ObsDn(rndI)] = cossgsea2(totalExpected, expectedDn, totalObserved, observedDn);

        end; 
        pvalueUp = sum(diffExp_ObsUp < realDifference)/length(diffExp_ObsUp);
        pvalueDn = sum(diffExp_ObsDn < realDifferenceDn)/length(diffExp_ObsDn);

        %disp(['Permutation test Up pValue: '  num2str(pvalueUp)]);
        %disp(['Permutation test Dn pValue: '  num2str(pvalueDn)]);
        %disp([' ']);
        fprintf(fileID,'\t%e\t%e\n',pvalueUp, pvalueDn);
      
        
   end; % for k
   
         for rndI = 1:1, % no for now...
        % Now a random miRNA
            myRandIdx = randperm(totalExpected);
            myRandIdx = myRandIdx(1:totalObserved);

            observedUp = 0;
            observedDn = 0;

            for i = 1:totalObserved,
                if (num(myRandIdx(1,i),indexFC-1) <= DnTh),
                    observedDn = observedDn +1;            
                    if (usePval ==1),
                        if (num(myRandIdx(1,i),PvalIndex) > pvalTh),
                            observedDn = observedDn -1;            
                        end;
                    end;
                end;
                if (num(myRandIdx(1,i),indexFC-1) > UpTh),
                    observedUp = observedUp +1;            
                    if (usePval ==1),
                        if (num(myRandIdx(1,i),PvalIndex) > pvalTh),
                            observedUp = observedUp -1;            
                        end;
                    end;
                end;
            end


            fprintf(fileID,'%s','Random Pathway');
            fprintf(fileID,'\t%f\t%f\t%f',totalExpected, expectedUp, expectedDn);
            fprintf(fileID,'\t%f\t%f\t%f',totalObserved, observedUp, observedDn);
            [pValUpReg, pValDnReg, pValDn_Up, ppWT, pmWT, pupWT, pdownWT, UpRegCILow, UpRegCIHigh, DnRegCILow, DnRegCIHigh] = cossgsea(totalExpected, expectedUp,expectedDn, totalObserved, observedUp, observedDn, 'Random', 0);
            fprintf(fileID,'\t%f\t%f\t%f\t%f\t%f\t%f\t%e\t%f\t%f\t%e',ppWT, pmWT, pupWT, pdownWT, UpRegCILow, UpRegCIHigh, pValUpReg, DnRegCILow, DnRegCIHigh, pValDnReg);
            fprintf(fileID,'\t%e\t%e\n',1, 1);

        end; 
       
        
   
   fclose(fileID);
   exit;
end; % for j
