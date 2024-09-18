function [F] = LoadFunSimpleFold(Node,U,icrm)
%% Compute current configuration 
Nodenw = Node;
Nodenw(:,1) = Node(:,1)+U(1:3:end);
Nodenw(:,2) = Node(:,2)+U(2:3:end);
Nodenw(:,3) = Node(:,3)+U(3:3:end);


LoadSize = 0.6;
%% Define Load here based on icrm and Nodenw
if icrm<=0
    error('Wrong increment!'); 
% elseif atan(Nodenw(4,1)/(Nodenw(4,3)+1e-10))<-0.1
%     DIRC = [Nodenw(4,3), 0, -Nodenw(4,1)];
%     Load = [4, LoadSize*DIRC/norm(DIRC)];
% elseif atan(Nodenw(4,1)/Nodenw(4,3))<0.1
%     DIRC = [1, 0, 0];
%     Load = [4, LoadSize*DIRC/norm(DIRC)];
else
    % DIRC3 = [-Nodenw(4,3), 0, -Nodenw(4,1)];
    DIRC4 = [Nodenw(4,3), 0, -Nodenw(4,1)];
    Load = [4, LoadSize*DIRC4/norm(DIRC4)];
end

%% Wrap up Load info to force/displacement vector F
m = size(Node,1);
F = zeros(3*m,1);
indp = Load(:,1);
F(3*indp-2) = Load(:,2);
F(3*indp-1) = Load(:,3);
F(3*indp) = Load(:,4);