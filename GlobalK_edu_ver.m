function [IF,K,state1,state2,state1_bend,state2_bend,state1_fold,state2_fold] = ...
    GlobalK_edu_ver(Ui,Node,truss,angles,state1,state2,state1_bend,state2_bend,state1_fold,state2_fold,iter,icrm)
Nn = size(Node,1);
IFb = zeros(3*Nn,1); IFp = IFb;
indi = zeros(36*size(truss.Bars,1),1); indj = indi; kentry = indi;
Nodenw(:,1) = Node(:,1)+Ui(1:3:end);
Nodenw(:,2) = Node(:,2)+Ui(2:3:end);
Nodenw(:,3) = Node(:,3)+Ui(3:3:end);

%state1_new=zeros(1,size(truss.Bars,1));
%state2_new=ones(1,size(truss.Bars,1));

% Bar
for bel = 1:size(truss.Bars,1) 
    eDof = [(-2:0)+(truss.Bars(bel,1)*3),(-2:0)+(truss.Bars(bel,2)*3)]';
    harden=state1(bel);
    pstr_p=state2(bel);
    % if (truss.Bars(bel, 1)>=7 && truss.Bars(bel, 2)>=7) || (truss.Bars(bel, 1)<=6 && truss.Bars(bel, 2)<=6)
    %     crease_flag = 0;
    % else
    %     crease_flag = 1;
    % end
    %===================================
    [~,Rbe,Kbe,harden,pstr_p] = BarKe(Ui(eDof),truss.B(bel,eDof),truss.L(bel),truss.A(bel),harden,pstr_p,truss.EMod(bel),truss.YBar(bel),truss.PMod(bel));
    % [~,Rbe,Kbe] = BarKe(Ui(eDof),truss.B(bel,eDof),truss.L(bel),truss.A(bel), iter);
     %[~,Rbe,Kbe,harden,pstr_p] = BarKe(Ui(eDof),truss.B(bel,eDof),truss.L(bel),truss.CM,truss.A(bel),harden,pstr_p);
    state1(bel)= harden;
    state2(bel)= pstr_p;
    %===================================
    IFb(eDof) = IFb(eDof)+Rbe;
    I=repmat(eDof,1,6); J=I';
    indi(36*(bel-1)+1:36*bel) = I(:);
    indj(36*(bel-1)+1:36*bel) = J(:); 
    kentry(36*(bel-1)+1:36*bel) = Kbe(:);
end
Kb = sparse(indi,indj,kentry,3*Nn,3*Nn);

% Hinge bend
indi = zeros(144*size(angles.bend,1),1); indj = indi; kentry = indi;
Lbend = truss.L(1:size(angles.bend,1));
for del = 1:size(angles.bend,1)
    eDof = reshape([3*angles.bend(del,:)-2;...
                    3*angles.bend(del,:)-1;...
                    3*angles.bend(del,:)],12,1);
    bend = angles.bend(del,:);
    harden_bend = state1_bend(del);
    pstr_p_bend = state2_bend(del);
    % [~,Rpe,Kpe] = FoldKe(Nodenw,bend,angles.Kb,angles.pb0(del),Lbend(del),angles.CMbend);
    [~,Rpe,Kpe,harden_bend,pstr_p_bend] = FoldKe_pl(Nodenw,bend,angles.Kb(del), ...
        angles.pb0(del),Lbend(del),harden_bend,pstr_p_bend,angles.ydb(del),angles.pmb(del),del,icrm);
    state1_bend(del)=harden_bend;
    state2_bend(del)=pstr_p_bend;
    IFp(eDof) = IFp(eDof)+Rpe;
    I=repmat(eDof,1,12); J=I';
    indi(144*(del-1)+1:144*del) = I(:);
    indj(144*(del-1)+1:144*del) = J(:); 
    kentry(144*(del-1)+1:144*del) = Kpe(:);
end;
Kbd = sparse(indi,indj,kentry,3*Nn,3*Nn);
if isempty(Kbd), Kbd = zeros(3*Nn); end

% Hinge fold
indi = zeros(144*size(angles.fold,1),1); indj = indi; kentry = indi;
Lfold = truss.L(size(angles.bend,1)+1:size(angles.bend,1)+size(angles.fold,1));
for fel = 1:size(angles.fold,1)
    eDof = reshape([3*angles.fold(fel,:)-2;...
                    3*angles.fold(fel,:)-1;...
                    3*angles.fold(fel,:)],12,1);
    fold = angles.fold(fel,:);
    harden_fold = state1_fold(fel);
    pstr_p_fold = state2_fold(fel);
    % [~,Rpe,Kpe] = FoldKe(Nodenw,fold,angles.Kf,angles.pf0(fel),Lfold(fel),angles.CMfold);
    [~,Rpe,Kpe,harden_fold,pstr_p_fold] = FoldKe_pl(Nodenw,fold,angles.Kf(fel), ...
        angles.pf0(fel),Lfold(fel),harden_fold,pstr_p_fold,angles.ydf(fel),angles.pmf(fel),fel,icrm);
    state1_fold(fel) = harden_fold;
    state2_fold(fel) = pstr_p_fold;
    IFp(eDof) = IFp(eDof)+Rpe;
    I=repmat(eDof,1,12); J=I';
    indi(144*(fel-1)+1:144*fel) = I(:);
    indj(144*(fel-1)+1:144*fel) = J(:); 
    kentry(144*(fel-1)+1:144*fel) = Kpe(:);
end;
Kfd = sparse(indi,indj,kentry,3*Nn,3*Nn);

IF = IFb+IFp;
K = Kb+Kbd+Kfd;
K = (K+K')/2;