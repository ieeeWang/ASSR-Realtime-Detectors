function result=PR2plot(x,y,plt,negtivelabel,ylim_range)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION
%   [auc]=AUC_pr(x,y)
% Plots the precision(PPV)-recall(sensitivity) curve and computes the area under the curve.
%
% INPUT ARGUMENTS:
%   x:      columun vector of data for both classes.
%   y:      column vector with the respective data labels. Each element
%           is equal to either -1 or 1.
%   plt:    if set to 1, a plot is generated. else value(e.g., 0): do not plot.
%   negtivelabel: default value is '-1', could be set to '0' according to y. 
% OUTPUT ARGUMENTS:
%   result.auc=auc; (area under precision-recall(PR) curve)
%   result.recall=recall;
%   result.precision=precision;
% (c) 2015 L. WANG @ TU/e
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
% Sort the rows in descending order using the values in column 1.
patterns = sortrows(patterns,-1);
y = patterns(:,2);
p = cumsum(y==1);
n = length(p);
% recall (or sensitivity)
recall = p/sum(y==1); 
% precision or positive predition value(PPV)
predict1=1:1:n; 
precision = p./predict1'; 

Y=(precision(2:n)+precision(1:n-1))/2;
X = recall(2:n) - recall(1:n-1);
auc=sum(Y.*X);
 
        
if (plt==1)
    figure
    plot(recall,precision,'-.k','LineWidth',2);
    grid on
    xlim([0 1])
    ylim(ylim_range)
    xlabel('Recall (sensitivity)','FontSize',16);ylabel('Precision (PPV)','FontSize',16);
    s=num2str(auc); s=strcat('AUC= ',s);
    title(['PR Curve, ' s],'FontSize',14);
end
    
result={};
result.auc=auc;
result.recall=recall;
result.precision=precision;
    