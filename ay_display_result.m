load('RES_NOV8_2_part_1.mat');
C = [COHAA];
F = [COHAA_F];

load('RES_NOV8_2_part_2.mat');
C = [C;COHAA];
F = [F;COHAA_F];

load('RES_NOV8_2_part_3.mat');
C = [C;COHAA];
F = [F;COHAA_F];

%% overall
load('RES_NOV8_2.mat');
TC = COHAA;
TF = COHAA_F;

%% Method
load('RES_NOV8_1.mat');
PC = COHA;
PF = COHA_F;

subplot(3,2,1)
imagesc(PC');
subplot(3,2,2)
imagesc(PF');

subplot(3,2,3)
imagesc(TC');
subplot(3,2,4)
imagesc(TF');


subplot(3,2,5)
imagesc(C');
subplot(3,2,6)
imagesc(F');
