function [MAT] = reGul(MAT,thr)
MAT(MAT>thr)=thr;
MAT(MAT<-thr)=-thr;
end