function [NODE,PANEL] = ConfigMorph(x_divs,y_divs,a,b,alfa1,alfa2,beta,sgngamma2)
if numel(alfa2) == 1, alfa2 = alfa2*ones(1,y_divs); end
if numel(sgngamma2) == 1, sgngamma2 = sgngamma2*ones(1,y_divs); end

alfa1 = alfa1/180*pi;
alfa2 = alfa2/180*pi;
beta = beta/180*pi;

cosgmma1 = cos(alfa1)/cos(beta);
cosgmma2 = cos(alfa2)/cos(beta);

numx = 2*x_divs; numy = 2*y_divs;
H1 = b*cosgmma1; H2 = a*cos(beta);
c = H1./cosgmma2;
H = H1+H2;
W = 2*a.*sin(beta);
L1 = b.*sqrt(1-cosgmma1.^2);
L2 = c.*sgngamma2.*sqrt(1-cosgmma2.^2);
L = L1+L2;

Ycoord = repmat([0;L1],1,y_divs)+repmat([0,cumsum(L(1:end-1))],2,1);
Ycoord = [Ycoord(:);sum(L)];
[X,Y] = meshgrid(linspace(0,W/2*numx,numx+1),Ycoord);

Z = 0*X;
Z(2:2:end,1:2:end) = H1; Z(1:2:end,2:2:end) = H2; 
Z(2:2:end,2:2:end) = H;
NODE = [reshape(X,[],1),reshape(Y,[],1),reshape(Z,[],1)];

k = 0; PANEL = cell(numx*numy,1);
for j=1:numy, for i=1:numx
        k = k+1;
        n1 = (i-1)*(numy+1)+j; n2 = i*(numy+1)+j;
        PANEL{k} = [n1 n2 n2+1 n1+1];
end, end