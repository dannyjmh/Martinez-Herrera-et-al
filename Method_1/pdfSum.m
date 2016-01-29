function [pdfZ,z]=pdfSum(pdfX,x,pdfY,y)
    delta=x(2)-x(1);
    [pdfZ,z]=myConv(pdfX,x,pdfY,y);
    pdfZ=pdfZ*delta;
end