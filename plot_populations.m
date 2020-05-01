clc %to clear the command
close all % to close eventual previous windows
clear all % to clear the workspace

FolderName = 'build/';

u1 = load([FolderName,'S.txt']);
u2 = load([FolderName,'Z.txt']);
u3 = load([FolderName,'R.txt']);
t = load([FolderName,'time.txt']);

plot(t, u1, 'b', t, u2, 'r', t, u3, 'g')
xlabel("Time");
ylabel("Individuals");

legend('S', 'Z', 'R')