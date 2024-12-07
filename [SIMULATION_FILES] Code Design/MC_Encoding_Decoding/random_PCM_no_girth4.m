function H = random_PCM_no_girth4(N,lambda,rho)

lambda=lambda(:)';
rho=rho(:)';
rand('twister',sum(1000000*clock));

degV_max = length(lambda);
degC_max = length(rho);

M = floor(sum(rho./(1:degC_max),2)/sum(lambda./(1:degV_max),2)*N);

preL = ((lambda./(1:degV_max)./sum(lambda./(1:degV_max),2))*N)';
cvx_begin
    cvx_quiet true 
%     cvx_solver sdpt3
%     cvx_solver SeDuMi
    cvx_solver Mosek
    variable L(degV_max) integer
    minimize( sum_square( L - preL ) )
    sum(L) == N
    L >= 0
cvx_end
L = (round(L))';

E = sum(L.*(1:degV_max),2); % E = (preN / sum(lambda./(1:degV_max),2));
preR = ((rho./(1:degC_max)./sum(rho./(1:degC_max),2))*M)';
degC_vec_inc = (1:degC_max);
cvx_begin
    cvx_quiet true 
%     cvx_solver sdpt3
%     cvx_solver SeDuMi
    cvx_solver Mosek
    variable R(degC_max) integer
    minimize( sum_square( R - preR ) )
    degC_vec_inc*R == E
    sum(R) == M
    R >= 0
cvx_end
R = (round(R))';


degV = zeros(1,N);
degC = zeros(1,M);


cnt=0;
for i=1:degV_max
    if L(1,i)>0
        for j=cnt+1:cnt+L(1,i)
            degV(1,j)=i;
        end
    cnt=j;
    end
end


cnt=0;
for i=1:degC_max
    if R(1,i)>0
        for j=cnt+1:cnt+R(1,i)
            degC(1,j)=i;
        end
    cnt=j;
    end
end

if sum(degV,2) ~= sum(degC,2)
    error('shit!');
end

% sock=0;
% for k = 1:N
%     for l = 1:Degv(k)
%         sock = sock +1;
%         sockvar(1,sock)=k;
%     end
% end

sock=0;
for k = 1:M
    for l = 1:degC(k)
        sock = sock +1;
        sockchk(1,sock)=k;
    end
end

% cnt = 0;
% for i = 1 : N
%     for j = 1 : Degv(i)
%         cnt = cnt +1;
%         VarSock(i,j)=cnt;
%     end
% end

cnt = 0;
for i = 1 : M
    for j = 1 : degC(i)
        cnt = cnt +1;
        ChkSock(i,j)=cnt;
    end
end

Nchk = zeros(M,degC_max);
Nvar = zeros(N,degV_max);
degvtemp = zeros(1,N);
degctemp = zeros(1,M);
selectsock = [];
H = false(M,N);

grDex = 0;
for i = N :-1:1
    i
    
    forbid_socket =[];
    candset = 1:sum(degV,2);
    candset = setdiff(candset,selectsock);
    for j = 1: degV(i)
        
        if j == 1
            %             rnd = randint(1,1,length(candset))+1;
            rnd = randi([0,length(candset)-1],1,1)+1;
        else
            %------------------Without Girth 4------------------%
            temp1 = Nchk(selectchk,:);
            Junkind1 = find(temp1==0);
            temp1(Junkind1)=[];
            
            temp3=[];
            for rr = 1 : length(temp1)
                temp2 = Nvar(temp1(rr),:);
                Junkind2 = find(temp2==0);
                temp2(Junkind2)=[];
                temp3 = [temp3,temp2];
            end
            temp4 = unique(temp3);
            temp5=[];
            for ss = 1 : length(temp4)
                tempsock = ChkSock(temp4(ss),:);
                Junkind3= find(tempsock==0);
                tempsock(Junkind3)=[];
                temp5=[temp5,tempsock];
            end
            forbid_socket = unique(temp5);
            candset = setdiff(candset,forbid_socket);
            %---------------------------------------------------------%
            %             rnd = randint(1,1,length(candset))+1;
            if isempty(candset)
                grDex = grDex+1;
                warning(grDex+" girth of length 4")
                candset = randi([1,size(sockchk,2)],1,1);
            end
            rnd = randi([0,length(candset)-1],1,1)+1;
        end
        selectsock = union(selectsock,candset(1,rnd));
        degvtemp(1,i)= degvtemp(1,i)+1;
        
        selectchk = sockchk(1,candset(1,rnd));
        degctemp(1,selectchk)=degctemp(1,selectchk)+1;
        
        H (selectchk , i) = 1;
        forbid_temp = ChkSock(selectchk,:);
        indzero = find(forbid_temp==0);
        forbid_temp(indzero)=[];
        forbid_socket=union(forbid_socket,forbid_temp);
        candset = setdiff(candset,forbid_socket);
        Nvar(i,j)=selectchk;
        Nchk(selectchk, degctemp(1,selectchk))=i;
        
    end
    forbid_socket = [];
    indzero = [];
    candset = [];
end