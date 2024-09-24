function [F] = LoadFunKresling(Node,U,icrm)
%% Compute current configuration 
Nodenw = Node;
Nodenw(:,1) = Node(:,1)+U(1:3:end);
Nodenw(:,2) = Node(:,2)+U(2:3:end);
Nodenw(:,3) = Node(:,3)+U(3:3:end);
DIRC = [0, 0, 1];
LoadSize = -0.003;
%% Define Load here based on icrm and Nodenw
if icrm<=0 
    error('Wrong increment!'); 
elseif icrm<=100
    Load = [7, LoadSize*DIRC/norm(DIRC);
            8, LoadSize*DIRC/norm(DIRC);
            9, LoadSize*DIRC/norm(DIRC);
            10, LoadSize*DIRC/norm(DIRC);
            11, LoadSize*DIRC/norm(DIRC);
            12, LoadSize*DIRC/norm(DIRC);];
elseif icrm<=250
    Load = [7, -LoadSize*DIRC/norm(DIRC);
            8, -LoadSize*DIRC/norm(DIRC);
            9, -LoadSize*DIRC/norm(DIRC);
            10, -LoadSize*DIRC/norm(DIRC);
            11, -LoadSize*DIRC/norm(DIRC);
            12, -LoadSize*DIRC/norm(DIRC);];
elseif icrm<=400
    DIRC = [0, 0, 1];
    Load = [7, LoadSize*DIRC/norm(DIRC);
            8, LoadSize*DIRC/norm(DIRC);
            9, LoadSize*DIRC/norm(DIRC);
            10, LoadSize*DIRC/norm(DIRC);
            11, LoadSize*DIRC/norm(DIRC);
            12, LoadSize*DIRC/norm(DIRC);];
elseif icrm<=550
    Load = [7, -LoadSize*DIRC/norm(DIRC);
            8, -LoadSize*DIRC/norm(DIRC);
            9, -LoadSize*DIRC/norm(DIRC);
            10, -LoadSize*DIRC/norm(DIRC);
            11, -LoadSize*DIRC/norm(DIRC);
            12, -LoadSize*DIRC/norm(DIRC);];
elseif icrm<=700
    DIRC = [0, 0, 1];
    Load = [7, LoadSize*DIRC/norm(DIRC);
            8, LoadSize*DIRC/norm(DIRC);
            9, LoadSize*DIRC/norm(DIRC);
            10, LoadSize*DIRC/norm(DIRC);
            11, LoadSize*DIRC/norm(DIRC);
            12, LoadSize*DIRC/norm(DIRC);];
elseif icrm<=850
    Load = [7, -LoadSize*DIRC/norm(DIRC);
            8, -LoadSize*DIRC/norm(DIRC);
            9, -LoadSize*DIRC/norm(DIRC);
            10, -LoadSize*DIRC/norm(DIRC);
            11, -LoadSize*DIRC/norm(DIRC);
            12, -LoadSize*DIRC/norm(DIRC);];
elseif icrm<=1000
    DIRC = [0, 0, 1];
    Load = [7, LoadSize*DIRC/norm(DIRC);
            8, LoadSize*DIRC/norm(DIRC);
            9, LoadSize*DIRC/norm(DIRC);
            10, LoadSize*DIRC/norm(DIRC);
            11, LoadSize*DIRC/norm(DIRC);
            12, LoadSize*DIRC/norm(DIRC);];
elseif icrm<=1450
    Load = [7, -LoadSize*DIRC/norm(DIRC);
            8, -LoadSize*DIRC/norm(DIRC);
            9, -LoadSize*DIRC/norm(DIRC);
            10, -LoadSize*DIRC/norm(DIRC);
            11, -LoadSize*DIRC/norm(DIRC);
            12, -LoadSize*DIRC/norm(DIRC);];
elseif icrm<=1650
    DIRC = [0, 0, 1];
    Load = [7, LoadSize*DIRC/norm(DIRC);
            8, LoadSize*DIRC/norm(DIRC);
            9, LoadSize*DIRC/norm(DIRC);
            10, LoadSize*DIRC/norm(DIRC);
            11, LoadSize*DIRC/norm(DIRC);
            12, LoadSize*DIRC/norm(DIRC);];
elseif icrm<=1850
    Load = [7, -LoadSize*DIRC/norm(DIRC);
            8, -LoadSize*DIRC/norm(DIRC);
            9, -LoadSize*DIRC/norm(DIRC);
            10, -LoadSize*DIRC/norm(DIRC);
            11, -LoadSize*DIRC/norm(DIRC);
            12, -LoadSize*DIRC/norm(DIRC);];
else
    DIRC = [0, 0, 1];
    Load = [7, LoadSize*DIRC/norm(DIRC);
            8, LoadSize*DIRC/norm(DIRC);
            9, LoadSize*DIRC/norm(DIRC);
            10, LoadSize*DIRC/norm(DIRC);
            11, LoadSize*DIRC/norm(DIRC);
            12, LoadSize*DIRC/norm(DIRC);];
end

%% Wrap up Load info to force/displacement vector F
m = size(Node,1);
F = zeros(3*m,1);
indp = Load(:,1);
F(3*indp-2) = Load(:,2); 
F(3*indp-1) = Load(:,3); 
F(3*indp) = Load(:,4);
