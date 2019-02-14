clear all;
close all;
clc;

javaclasspath('../RCA1Portable/jars/tda.jar');
import api.*;
tda=Tda();


N = 50;
RNG = 12325;
rng(RNG);
[X Y] = gen_annulus(1,2,N);

plot(X,Y,'.');
axis equal;

LIM = 2+1.25;
RS = [0.1 0.32 0.38 .65 1.25];
figure;
CM = zeros(4);
load('cm');
for i=1:length(RS)
    figure;
    R = RS(i);
    disks(X,Y,R,100, 0.8*[1 1 1]);
    hold on;
%     scatter(X,Y,50,'k');
    plot(X,Y,'ok', 'linewidth',0.5,'markerfacecolor','k' )
    xlim([-LIM LIM]);
    ylim([-LIM LIM]);
    axis square;
    axis off;
    sfname = sprintf('Output/ph_example_%d',i);
   %print(sfname, '-dpng');
end

MN = 0;
MX = 1.5;
TH = 0.01;
COLORS  ={[191 0 191]/256, [204 255 51]/256} ;

%%
tda.RCA1( { 'settingsFile=../RCA1Portable/data/cts.txt', 'supplyDataAs=pointCloud', 'distanceBoundOnEdges=5'}, [X'  Y']);
D = {};
B = {};
c = 1;
J = tda.getResultsRCA1(0).getIntervals;
BS{1} = J(:,1);
tmp = J(:,2);
idx = find(tmp == -1);
tmp(idx) = Inf;
DS{1} = tmp;
J = tda.getResultsRCA1(1).getIntervals;
if (length(J) == 0)
    B{2} = [];
    D{2} = [];
    return;
end
BS{2} = J(:,1);
tmp = J(:,2);
idx = find(tmp == -1);
tmp(idx) = Inf;
DS{2} = tmp;
%%
ND = 2;

figure;

for i=1:ND
    B = BS{i}/2;
    D = DS{i}/2;
    
    idx_inf = find(D == inf);
    D(idx_inf) = MX;
    
    if (i == 1)
        [D idx] = sort(D,'descend');
        B = B(idx);
    else
        [B idx] = sort(B);
        D = D(idx);
    end
    
    for j=1:length(B)
        if (D(j)-B(j) > TH)
            plot([B(j) D(j)], [c c], '-', 'Color', COLORS{i},  'linewidth', 5);
            hold on;
            c = c+1;
        end
    end
    if (i < ND)
        plot([MN MX], [c c], 'k--');
        c = c+1;
    end
end
set(gca, 'fontsize', 20);
set(gca, 'ydir', 'reverse');
set(gca, 'ytick', []);
ylim([-0.5 c+0.5]);
xlim([MN MX]);
hold off;
%print('Output/ph_example_barcode', '-dpng');

%%
% D = DS{1};
% B = BS{1};
% mx = max(D(find(D< inf)));
% idx = find(D-B > TH & D < inf);
% idx_inf = find(D == inf);
% figure;
% plot(B(idx),D(idx),'.');
% hold on;
D = DS{2}/2;
B = BS{2}/2;
% mx = max(mx, max(D));
mx = max(D);
idx = find(D-B > TH);
figure;
plot(B(idx),D(idx),'ok', 'markersize', 15,'linewidth', 3,'markerfacecolor', CM(32,:));
hold on;
mx = 1.2*mx;
plot([0 mx], [0 mx],'k');
axis equal
xlim([0 mx]);
ylim([0 mx]);
xlabel('birth');
ylabel('death');
set(gca, 'fontsize', 20);
set(gca, 'ytick', 0:0.2:mx);
set(gca, 'xtick', 0:0.2:mx);
%print('Output/ph_example_diagram', '-dpng');


