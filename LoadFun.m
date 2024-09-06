function [F] = LoadFun(Node,U,icrm)
%% Compute current configuration 
Nodenw = Node;
Nodenw(:,1) = Node(:,1)+U(1:3:end);
Nodenw(:,2) = Node(:,2)+U(2:3:end);
Nodenw(:,3) = Node(:,3)+U(3:3:end);
DIRC = [0, 1, 0];
Mult = 0.0002;
PieceStep = 200;
%% Define Load here based on icrm and Nodenw
if icrm<=0 
    error('Wrong increment!'); 
elseif icrm<=PieceStep
    Load = [3, Mult*DIRC/norm(DIRC)];
elseif icrm<=3*PieceStep
    Load = [3, -Mult*DIRC/norm(DIRC)];
elseif icrm<=5*PieceStep
    Load = [3, Mult*DIRC/norm(DIRC)];
elseif icrm<=7*PieceStep
    Load = [3, -Mult*DIRC/norm(DIRC)];
elseif icrm<=9*PieceStep
    Load = [3, Mult*DIRC/norm(DIRC)];
elseif icrm<=11*PieceStep
    Load = [3, -Mult*DIRC/norm(DIRC)];
else
    Load = [3, Mult*DIRC/norm(DIRC)];
end

%% Wrap up Load info to force/displacement vector F
m = size(Node,1);
F = zeros(3*m,1);
indp = Load(:,1);
F(3*indp-2) = Load(:,2); 
F(3*indp-1) = Load(:,3); 
F(3*indp) = Load(:,4);
