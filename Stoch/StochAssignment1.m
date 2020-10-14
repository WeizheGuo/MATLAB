%Stoch Assignment 1
%Weizhe Guo
%DND

%% 1)
% a)
trials = 100000;

sim_3d6 = sum(randi(6,[3,trials]));
prob_18 = nnz(sim_3d6 == 18)/trials

% b)
sim_3d6_fun = max(sum(randi(6,[3,3,trials])),[],2);
% take max along the second dimension to get the maximum of three rolls
prob_fun_18 = nnz(sim_3d6_fun(1,1,:)==18)/trials

% c)
sim_3d6_sixabi =squeeze(max(sum(randi(6,[3,3,6,trials])),[],2));
prob_sixabimax = nnz(~any(sim_3d6_sixabi-18,1))/trials
% Thanks for Alex Zheng Liu helping on this
% This method avoids doing for loops and is thus very efficient

% d)
prob_sixabiavg = nnz(~any(sim_3d6_sixabi-9,1))/trials

%% 2)
% a)
sim_1d4 = randi(4,[1,trials]);
avg_hp = mean(sim_1d4)
sim_2d2 = sum(randi(2,[2,trials]));
avg_dam = mean(sim_2d2)
prob_3dam = nnz(sim_2d2>3)/trials

% b)
prob_2dam = nnz(sim_2d2==2)/trials
prob_3dam = nnz(sim_2d2==3)/trials
prob_4dam = nnz(sim_2d2==4)/trials

prob_1hp = nnz(sim_1d4==1)/trials
prob_2hp = nnz(sim_1d4==2)/trials
prob_3hp = nnz(sim_1d4==3)/trials
prob_4hp = nnz(sim_1d4==4)/trials

% c)
sim_six_trolls = randi(4,[6,trials]);
prob_destroyall = nnz(all(sim_six_trolls <= sim_2d2,1))/trials

% d)
beat_num = sum(sim_six_trolls<=sim_2d2,1);
max_sim_six_trolls = max(sim_six_trolls,[],1);
num_occur = 0;
total_hp = 0;

for i = 1:trials
    if beat_num(1,i) == 5
        total_hp = total_hp + max_sim_six_trolls(1,i);
        num_occur = num_occur + 1;
    end
end

avg_hp = total_hp/num_occur
% I calculated the hp before the fire ball hits the troll.

% e)
sword_dam = sum(randi(6,[2,trials]));
hammer_dam = randi(4,[1,trials]);
total_dam = 0;

attack_condition = randi(2,[2,trials]);
for i = 1:trials
    if attack_condition(1,i)==2
        total_dam = total_dam + sword_dam(1,i);
        if attack_condition(2,i)==2
        total_dam = total_dam + hammer_dam(1,i);
        end
    end
end

avg_dam = total_dam/trials

% f)
% I believe you can! So the prob is 1.