clc; clearvars; close all;
addpath('C:\Users\Bart\Documents\5.4 integration project 5SC26\GitHubFolderStuff\5sc26_group7\Model');
%%
c46_x1_optimum = [5.53388313158027e-05,0.628942248124698,0.828797041915195,...
    0.942461739693451,0.932252052511358,1.02442216498716,0.995209299898076,...
    1.01699750453345,1.05393097018382,0.977373328855299,0.955416378269874];
c46_x2_optimum = [6.20619791555931e-05,0.918601609759064,0.847702416510482,...
    1.10488850638160,1.07685278368839,1.10549728448027,1.12223314407549,...
    1.12778732487323,1.15249953736365,1.09080730625011,1.05869028025712];

u_sampled = 0.01:0.01:0.1;
u1_poly = Pump1Function(u_sampled);
u2_poly = Pump2Function(u_sampled);


figure(1)
set(gca,'LooseInset',get(gca,'TightInset'));
hold on
a=plot(u_sampled,c46_x1_optimum(2:11).*u_sampled,'linestyle','none','marker',...
    'o','MarkerFaceColor','green','MarkerEdgeColor','green','Markersize',6)
b=plot(u_sampled,c46_x2_optimum(2:11).*u_sampled,'linestyle','none','marker',...
    'o','MarkerFaceColor','blue','MarkerEdgeColor','blue','Markersize',6)
c=plot(0:0.01:0.1,0:0.01:0.1,'linestyle','--','color','black')
d=plot(0.01:0.01:0.1,u1_poly,'linestyle','--','color','green')
e=plot(0.01:0.01:0.1,u2_poly,'linestyle','--','color','blue')
legend([a b c d e],'Fmincon fit pump 1','Fmincon fit pump 2','expected value','polyfit pump 1','polyfit pump 2','location','northwest')
xlim([0 0.1])
xlabel('Q_{theoretic} [L/s]')
ylim([0 0.11])
ylabel('Q_{actual} [L/s]')

grid on