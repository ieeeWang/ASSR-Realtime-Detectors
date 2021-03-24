function result=ROC2plot(x,y,plt,negtivelabel,ylim_range)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION
%   [auc]=ROC(x,y)
% Plots the ROC curve and computes the area under the curve.
%
% INPUT ARGUMENTS:
%   x:      column vector of data for both classes.
%   y:      column vector with the respective data labels. Each element
%           is equal to either -1 or 1.
%   plt:    if set to 1, a plot is generated.
% OUTPUT ARGUMENTS:
%   result.auc:    area under ROC curve.
%   result.recall=recall;
%   result.precision=precision;
% (c) 2010 S. Theodoridis, A. Pikrakis, K. Koutroumbas, D. Cavouras
% update by Lei@TUe 2016 Mar. 22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin>=4
    negtivevalue = negtivelabel;
else
    negtivevalue = -1; % default value
end

% ensure it to columun
x = x(:);
y = y(:);
patterns = [x y];
patterns = sortrows(patterns,-1); % Sort the rows in descending order using the values in column 1.
y = patterns(:,2);
p = cumsum(y==1);
tp = p/sum(y==1); % true positive/all positive -> sensitivity
n = cumsum(y==negtivevalue);
fp = n/sum(y==negtivevalue);% false negtive/all negtive -> (1- specificity)

n = length(tp);
Y=(tp(2:n)+tp(1:n-1))/2;
X = fp(2:n) - fp(1:n-1);
% auc=sum(Y.*X)-0.5;
auc=sum(Y.*X);

if (plt==1)
    plot(fp,tp,'-.k','LineWidth',2);
    grid on
    xlim([0 1])
    ylim(ylim_range)
    xlabel('(1-specificity)','FontSize',16);ylabel('Sensitivity','FontSize',16);
    s=num2str(auc); s=strcat('AUC= ',s);
    title(['ROC Curve, ' s],'FontSize',14);hold off;
end

result={};
result.auc=auc;
result.sensitivity=tp;
result.FPrate=fp;% False Positive rate = (1- specificity)
