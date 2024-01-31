function [x_roc, y_roc, y_prc, auc_roc,auc_prc, cutoff_opt_ba, balanced_accuracy_opt_ba,...
    sensitivity_opt_ba, specificity_opt_ba,precision_opt_ba,accuracy_opt_ba] = performance_wBA(pred,Y)

%performance metrics with balanced accuracy (BA) included
% Input:
% - pred: vector of predicited probabilities
% - Y: vector of the outcome
% Output:
% - x_roc: roc curve x axis
% - y_roc: roc curve y axis
% - y_prc: y_axis for precision-recall curve
% - auc_roc: area under the roc curve
% - cutoff_opt_roc: cutoff that maximizes accuracy
% - cutoff_opt_prc: cutoff that minimizes difference between specificity and
%   sensitivity
% - accuracy_opt_roc: accuracy at cutoff_opt_roc
% - accuracy_opt_prc: accuracy at cutoff_opt_prc
% - sensitivity_opt_roc: sensitivity at cutoff_opt_roc
% - sensitivity_opt_prc: sensitivity at cutoff_opt_prc
% - specificity_opt_roc: specificity at cutoff_opt_roc
% - specificity_opt_prc: specificity at cutoff_opt_prc
% - precision_opt_roc: precision at at cutoff_opt_roc
% - precision_opt_prc: precision at at cutoff_opt_prc
% - auc_prc: area under the precision-recall curve

cutoff_vec = sort(pred,'descend');
sensitivity = zeros(size(cutoff_vec));
specificity = zeros(size(cutoff_vec));
accuracy = zeros(size(cutoff_vec));
balanced_accuracy = zeros(size(cutoff_vec)); 
precision=zeros(size(cutoff_vec));


for k = 1:length(cutoff_vec)
    
    c= cutoff_vec(k);
    
    TP = length(find(pred>=c & Y==1));
    FP = length(find(pred>=c & Y==0));
    TN = length(find(pred<c & Y==0));
    FN = length(find(pred<c & Y==1));
    
    sensitivity(k) = TP/(TP+FN);
    specificity(k) = TN/(TN+FP);
    accuracy(k) = (TP+TN)/(TP+FP+TN+FN);
    precision(k)=TP/(TP+FP);
    balanced_accuracy(k)=0.5*((TP/(TP+FN))+(TN/(TN+FP)));

end


x_roc = 1-specificity;
y_roc = sensitivity;
y_prc = precision;

auc_roc = trapz(x_roc, y_roc);
auc_prc = trapz(y_roc, y_prc);

%Cutoff for maximizing balanced accuracy
Iopt_ba= find(balanced_accuracy==max(balanced_accuracy));
cutoff_opt_ba = cutoff_vec(Iopt_ba);

sensitivity_opt_ba = sensitivity(Iopt_ba);
specificity_opt_ba = specificity(Iopt_ba);
balanced_accuracy_opt_ba= balanced_accuracy(Iopt_ba);
accuracy_opt_ba=accuracy(Iopt_ba);
precision_opt_ba=precision(Iopt_ba);

% Cutoff for maximizing accuracy
Iopt_roc = find(accuracy==max(accuracy),1,'last');
cutoff_opt_roc = cutoff_vec(Iopt_roc);


% Cutoff for minimizing difference between specificity and sensitivity
Iopt_prc = find(abs(specificity-sensitivity)==min(abs(specificity-sensitivity)),1,'last');
cutoff_opt_prc=cutoff_vec(Iopt_prc);



end