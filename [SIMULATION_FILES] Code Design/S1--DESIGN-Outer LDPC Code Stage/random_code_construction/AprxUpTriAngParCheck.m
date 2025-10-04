function [H_UpTriAng,g,VecCol,VecRow]=AprxUpTriAngParCheck(H)

% This function simulate Approximate Upper Triangular Parity-Check Matrix
% Proposed By T. Richardson & R. Urbanke 
% This Program Written By: R. Asvadi PhD Candidate K. N. Toosi Univ of Tech
% in 9/22/2008 ver: 0.1

[M,N]=size(H);
% rnk=gfrank(H);
% if rnk< M
%     error('Input Matrix Must be Full Row Rank to Program Work Correctly')
% end

K=N-M;
%Initialization
t=0;
g=0;
% Htmp=H;
% clear H;
VecRow=1:M;%This Vector Save Original index of H rows.
VecCol=1:N;%This Vector Save Original index of H columns.
while(t< N-K-g)
    [t N-K-g]
    cnt=0;
    ResDeg=zeros(1,N-t);
    %     for j=t+1:c
    %         Sum=0;
    %         for i=t+1:N-K-g
    %             Sum=Sum+Htmp(i,j);
    %         end
    %         ResDeg(j-t,1)=Sum;
    %     end
    Hresedual=H(t+1:N-K-g,t+1:N);
    ResDeg=sum(Hresedual,1);
    [minval,minind]=min(ResDeg);
    ResDegTmp=ResDeg;
    while (minval==0)&& (cnt<N)
        cnt=cnt+1;
        ResDegTmp(1,minind)=999;
        [minval,minind]=min(ResDegTmp);
    end
    MinResDeg=minval;
    MinResDegInd=find(ResDeg==minval);

    if (MinResDeg==1)%Extend in Algorithm
        rndint=randi(length(MinResDegInd),1,1);
%         rndint=1;
        col=t+MinResDegInd(1,rndint);
        %         cntrow=0;
        %         for r=t+1:N-K-g
        %             if (Htmp(r,col)~=0)
        %                 cntrow=cntrow+1;
        %                 row(cntrow,1)=r;
        %             end
        %         end
        row=t+find(Hresedual(:,MinResDegInd(1,rndint))==1);

        tmpvar0=VecCol(1,col);
        tmpvar1=VecCol(1,t+1);
        VecCol(1,t+1)=tmpvar0;
        VecCol(1,col)=tmpvar1;

        vectmpc=H(:,col);%Swap Column col with column t+1
        H(:,col)=H(:,t+1);
        H(:,t+1)=vectmpc;

        if isempty(row)==1
            error('There has been an error accured')
        end

        vectmpr=H(row(1,1),:);%Swap row r1 with row t+1
        H(row(1,1),:)=H(t+1,:);
        H(t+1,:)=vectmpr;

        tmpvar0=VecRow(1,row(1,1));
        tmpvar1=VecRow(1,t+1);
        VecRow(1,t+1)=tmpvar0;
        VecRow(1,row(1,1))=tmpvar1;


        t=t+1;

    else %Choose in Algorithm
        rndint=randi(length(MinResDegInd),1,1);
%         rndint=1;
        col=t+MinResDegInd(1,rndint);%Associate to original column index!!
        %         cntrow=0;
        %         for r=t+1:N-K-g
        %             if (Htmp(r,col)~=0)
        %                 cntrow=cntrow+1;
        %                 row(cntrow,1)=r;
        %             end
        %         end

        %         for r=1:N-K-g-t
        %             if (Hresedual(r,MinResDegInd(rndint,1))~=0)
        %                 cntrow=cntrow+1;
        %                 row(cntrow,1)=;
        %             end
        %         end
        row=t+find(Hresedual(:,MinResDegInd(1,rndint))==1);


        tmpvar0=VecCol(1,col);
        tmpvar1=VecCol(1,t+1);
        VecCol(1,t+1)=tmpvar0;
        VecCol(1,col)=tmpvar1;

        vectmpc=H(:,col);%Swap Column col with column t+1
        H(:,col)=H(:,t+1);
        H(:,t+1)=vectmpc;

        vectmpr=H(row(1,1),:);%Swap row r1 with row t+1
        H(row(1,1),:)=H(t+1,:);
        H(t+1,:)=vectmpr;

        tmpvar0=VecRow(1,row(1,1));
        tmpvar1=VecRow(1,t+1);
        VecRow(1,t+1)=tmpvar0;
        VecRow(1,row(1,1))=tmpvar1;

        for k=MinResDeg:-1:2  %move rows r2,r3,...rd to the bottom of the matrix
            vectmpr=H(row(k,1),:);%Move rd,...,r2 to bottom of the matrix
            H(row(k,1),:)=H(M-g+k-MinResDeg,:);
            H(M-g+k-MinResDeg,:)=vectmpr;

            tmpvar0=VecRow(1,row(k,1));
            tmpvar1=VecRow(1,M-g+k-MinResDeg);
            VecRow(1,row(k,1))=tmpvar1;
            VecRow(1,M-g+k-MinResDeg)=tmpvar0;

        end
        t=t+1;
        g=g+MinResDeg-1;

    end
    clear ResDeg ResDegTmp MinResDeg MinResDegInd vectmpc vectmpr Hresedual row col
end
%------------------Encoder----------------------------%
H_UpTriAng = H;