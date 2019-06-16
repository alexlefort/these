clear all
close all
clc

p = ureal('p',1,'range',[0 10]);
k = realp('k',1);
k.Minimum = 0;
k.Maximum = 10;

sys = tf(1, [1/(p+1) (1+k)]);

disp(sys.Blocks.p);
sys.Blocks.p = ureal('p',2,'range',[1 10]);
disp(sys.Blocks.p);
set(p,'range',[0.5 20]);
disp(p);
disp(sys.Blocks.p);