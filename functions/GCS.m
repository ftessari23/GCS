function [GCS_mean,GCS_vec] = GCS(d_avgp,d_maxp,d_minp,plot_flag)
%Geometric Configuration Similarity (GCS)

%Maximum Tip Fingers Distances ordered as: (THB, IDX, MID, RIG, LIT) x (TIP,DIP,PIP)
load('mean_d_max_15_final.mat');

%Compute GCS
for i = 1:15
    for j = 1:3
        if j == 1
            GCS_vec(j,i) = 100*(1-d_avgp(i)/mean_d_max(i));
        elseif j == 2
            GCS_vec(j,i) = 100*(1-d_maxp(i)/mean_d_max(i));
        elseif j == 3
            GCS_vec(j,i) = 100*(1-d_minp(i)/mean_d_max(i));
        end
    end
end

%Average for TIP, DIP, PIP
for j = 1:3
    GCS_vec_favg(j,1) = mean([GCS_vec(j,1),GCS_vec(j,6),GCS_vec(j,11)]);
    GCS_vec_favg(j,2) = mean([GCS_vec(j,2),GCS_vec(j,7),GCS_vec(j,12)]);
    GCS_vec_favg(j,3) = mean([GCS_vec(j,3),GCS_vec(j,8),GCS_vec(j,13)]);
    GCS_vec_favg(j,4) = mean([GCS_vec(j,4),GCS_vec(j,9),GCS_vec(j,14)]);
    GCS_vec_favg(j,5) = mean([GCS_vec(j,5),GCS_vec(j,10),GCS_vec(j,15)]);

    GCS_vec_fstd(j,1) = std([GCS_vec(j,1),GCS_vec(j,6),GCS_vec(j,11)]);
    GCS_vec_fstd(j,2) = std([GCS_vec(j,2),GCS_vec(j,7),GCS_vec(j,12)]);
    GCS_vec_fstd(j,3) = std([GCS_vec(j,3),GCS_vec(j,8),GCS_vec(j,13)]);
    GCS_vec_fstd(j,4) = std([GCS_vec(j,4),GCS_vec(j,9),GCS_vec(j,14)]);
    GCS_vec_fstd(j,5) = std([GCS_vec(j,5),GCS_vec(j,10),GCS_vec(j,15)]);
end

GCS_mean = mean(GCS_vec_favg,2);

%Plot colors
color_idx = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.4660 0.6740 0.1880];

if plot_flag == 1
    figure('color','white'),
    % tt = tiledlayout(1,3);
    for j = 1:3
        % nexttile(j)
        polarplot((pi/180)*[0 72 144 216 288 360],[GCS_vec_favg(j,:),GCS_vec_favg(j,1)],'o-','LineWidth',2,'MarkerSize',8,'Color',color_idx(j,:),'MarkerFaceColor',color_idx(j,:)),
        hold on
        thetaticks(0:72:360)
        rlim([0 100])
        thetaticklabels({'Thumb','Index','Middle','Ring','Little'})
        set(gca,'FontName','Times New Roman')
        % title(strcat('mean(GCS)=',num2str(round(GCS_mean(j),1)),'%'))
        % if j == 1
        %     subtitle('Average Posture')
        % elseif j == 2
        %     subtitle('Max Posture')
        % elseif j == 3
        %     subtitle('Min Posture')
        % end
    end
    legend(strcat('Avg. Posture, mean(GCS)=',num2str(round(GCS_mean(1),1)),'\%'),...
           strcat('Max. Posture, mean(GCS)=',num2str(round(GCS_mean(2),1)),'\%'),...
           strcat('Min. Posture, mean(GCS)=',num2str(round(GCS_mean(3),1)),'\%'),'Location','southoutside','FontName','Times New Roman');
    title('Geometric Configuration Similarity');
    % set(gca,'FontName','TimesNewRoman')
end
end