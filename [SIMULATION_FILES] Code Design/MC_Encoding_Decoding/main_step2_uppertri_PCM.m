clc
clearvars

load PCM_1.mat H; % load the parity check matrix

[M,N] = size(H); K = N-M;
[H_ut,g,permDex]=AprxUpTriAngParCheck(H);

T = H_ut(1:N-K-g,1:N-K-g);
A = H_ut(1:N-K-g,N-K-g+1:N-K);
B = H_ut(1:N-K-g,N-K+1:N);
E = H_ut(N-K-g+1:N-K,1:N-K-g);
C = H_ut(N-K-g+1:N-K,N-K-g+1:N-K);
D = H_ut(N-K-g+1:N-K,N-K+1:N);

gfrank(T,2) == size(T,1)
inv_T = inv_gfmat(T);

preMulmat = logical([eye(N-K-g),zeros(N-K-g,g);mod(-E*inv_T,2),eye(g)]);

encMtx = logical(mod(preMulmat*H_ut,2));
phi = encMtx(N-K-g+1:N-K,N-K-g+1:N-K);

if gfrank(phi,2) ~= size(phi,1)
    warning("The 'phi' matrix does not have a full rank, so this part tries to make it " + ...
        "full-rank by column permutations; It might take a while and might not succeed " + ...
        "at the end! If u don't receive further warning messages or errors, everything's " + ...
        "gonne be fine, don't freak out!")
    % input('Press ''Enter'' to continue...','s');
    clear B D encMtx phi
    C_perm = zeros(2,g);
    C_perm(1,:) = (N-K-g+1:N-K);

    auxC = false(size(C));
    auxA = false(size(A));

    [srtpotaddex_val,srtpotaddex] = sort( sum(H_ut(N-K-g+1:N-K,:),1) );

    potaddexDex = length(srtpotaddex);
    % auxC(:,1) = H_ut(N-K-g+1:N-K,srtpotaddex(potaddexDex));
    % auxA(:,1) = H_ut(1:N-K-g,srtpotaddex(potaddexDex));
    % C_perm(2,1) = srtpotaddex(potaddexDex);

    auxphi = mod(auxC+(preMulmat(N-K-g+1:end,1:N-K-g))*auxA,2);
    phiRank = gfrank(auxphi,2);

    cDex=1;
    % for srtpotaddex = length(srtpotaddex)-1:-1:max(find(srtpotaddex_val==0))+1
        while gfrank(auxphi,2) ~= size(auxphi,1)
            auxC(:,cDex)=H_ut(N-K-g+1:N-K,srtpotaddex(potaddexDex));
            auxA(:,cDex)=H_ut(1:N-K-g,srtpotaddex(potaddexDex));
            auxphi = mod(auxC+(preMulmat(N-K-g+1:end,1:N-K-g))*auxA,2);
            if gfrank(auxphi,2)>phiRank
                C_perm(2,cDex) = srtpotaddex(potaddexDex);
                phiRank = gfrank(auxphi,2);
                [phiRank g]
                cDex = cDex+1;
            end
            if cDex > g %|| potaddexDex < max(find(srtpotaddex_val==0))+1
                continue
            end
            potaddexDex=potaddexDex-1;
            if potaddexDex == 0
                warning("This (specific) 'H' matrix can not be used for encoding by this method; " + ...
                    "if you want to proceed, this program will remove " ...
                    + num2str(g-phiRank) + " rows/cols from 'H_ut' and " + ...
                    num2str(g-phiRank) + " cols from 'H', so the " + ...
                    "generated codewords might be noisy w.r.t. 'H'. " + ...
                    "(You can use the uppertriangularized matrix 'H_ut'" + ...
                    " along with systematic codewords 'codeword_sys' for decoding---" + ...
                    "'codeword_sys' is not noisy w.r.t. 'H_ut'. However, notice that 'H_ut' " + ...
                    "might not be sparse and obviously include a lot of short cycles.) " + ...
                    "BTW check the syndromes!")
                    % input('Press ''Enter'' to continue...','s');
                break
            end
        end
% end
    if gfrank(auxphi,2) == size(auxphi,1)
        C = auxC; clear auxC;
        A = auxA; clear auxA;
        H_ut(:,C_perm(2,:)) = H_ut(:,C_perm(1,:));
        H_ut(:,C_perm(1,:)) = [A;C];

        auxpermDex = permDex(C_perm(2,:));
        permDex(C_perm(2,:)) = permDex(C_perm(1,:));
        permDex(C_perm(1,:)) = auxpermDex;

        B = H_ut(1:N-K-g,N-K+1:N);
        D = H_ut(N-K-g+1:N-K,N-K+1:N);
        encMtx = logical(mod(preMulmat*H_ut,2));
        phi = encMtx(N-K-g+1:N-K,N-K-g+1:N-K);
    else
        newg = gfrank(auxphi,2);
        auxA(:,newg+1:end)=[];
        auxC(:,newg+1:end)=[];
        H(:,permDex(N-K-g+newg+1:N-K))=[];
        permDex(N-K-g+newg+1:N-K)=[];

        newE = false(newg, size(E,2));
        newC = false(newg);

        D = H_ut(N-K-g+1:N-K,N-K+1:N);
        newD = false(newg, size(D,2));

        newphi = mod(newC-newE*inv_T*auxA,2);
        newphiRank = gfrank(newphi,2);

        cnt = 1;
        for oldphi_rDex = 1:g
            newC(cnt,:)=auxC(oldphi_rDex,:);
            newE(cnt,:)=   E(oldphi_rDex,:);
            newD(cnt,:)=   D(oldphi_rDex,:);

            newphi = mod(newC-newE*inv_T*auxA,2);

            if gfrank(newphi,2)>newphiRank
                [newphiRank newg]
                newphiRank = gfrank(newphi,2);
                cnt=cnt+1;
            end

            if cnt>newg
                continue
            end
        end
        newB = H_ut(1:N-K-g,N-K+1:N); clear H_ut;
        H_ut = [   T, auxA, newB;...
                newE, newC, newD];
        clear A B E C D g
        [M,N] = size(H_ut); K = N-M;
        g = newg;

        A = H_ut(1:N-K-g,N-K-g+1:N-K);
        B = H_ut(1:N-K-g,N-K+1:N);
        E = H_ut(N-K-g+1:N-K,1:N-K-g);
        C = H_ut(N-K-g+1:N-K,N-K-g+1:N-K);
        D = H_ut(N-K-g+1:N-K,N-K+1:N);

        preMulmat = logical([eye(N-K-g),zeros(N-K-g,g);mod(-E*inv_T,2),eye(g)]);

        encMtx = logical(mod(preMulmat*H_ut,2));
        phi = encMtx(N-K-g+1:N-K,N-K-g+1:N-K);

    end

end

inv_phi = inv_gfmat(phi);

preCalc_p2=logical(mod(mod(-inv_phi,2)*mod((D-E*inv_T*B),2),2));

%%

s = (rand(1,K) < 0.5);
p2=logical(mod(preCalc_p2*s',2));
p1=logical(mod(-inv_T*(mod(A*p2,2)+mod(B*s',2)),2));
codeword_sys=[p1;p2;s'];


[truDex_val,truDex]=sort(permDex);
codeword(1:N)=codeword_sys(truDex);
codeword=codeword';

sum(mod(H*codeword,2))

save ENC.mat -v7.3