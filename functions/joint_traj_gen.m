function joint_trj = joint_traj_gen(w,w_avg,tf,ts,minA,maxA,opt)

%Fixed Parameters
if opt == 1
    %For Nice Graphix
    t = 0:ts:tf;
    joint_trj(:,1) = t;
    p = polyfit([0 tf/4 tf/2 3*tf/4 tf],[0 maxA 0 minA 0],4);
    u = p(1)*t.^4+p(2)*t.^3+p(3)*t.^2+p(4)*t+p(5);
else
    %For Computation
    t = [0 ts 2*ts 3*ts 4*ts];
    joint_trj(:,1) = t;
    u = [0 maxA minA 0 0];
end

for i = 1:length(w)
    if (i == 4) || (i == 8) || (i == 12) || (i == 16) || (i == 20)
        joint_trj(:,i+1) = +w(i)*u+w_avg(i);
    else
        joint_trj(:,i+1) = +w(i)*u+w_avg(i);
    end
end

end