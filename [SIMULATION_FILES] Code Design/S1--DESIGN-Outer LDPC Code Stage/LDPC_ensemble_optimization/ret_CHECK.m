clear;clc;
n_acc=inf;
n=inf;
for iter=0:113
    load(num2str(iter));
    if n<n_acc
        iter_acc=iter;
        n_acc=n;
    end
end