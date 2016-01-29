function [pdfZ,z]=pdfDiff(pdfX,x,pdfY,y)
    delta=x(2)-x(1);
    [yp,idx]=sort(-y);
    pdfYp=pdfY(idx);
    [pdfZ,z]=myConv(pdfX,x,pdfYp,yp);
    pdfZ=pdfZ*delta;
end