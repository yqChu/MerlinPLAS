function stopflag = StoppingCriteria(Node,U,icrm,Fhis)
% stopflag = (-U(121*3-2))>12;
a = 1;
stopflag = icrm>100 && Fhis(1)<1e-3;