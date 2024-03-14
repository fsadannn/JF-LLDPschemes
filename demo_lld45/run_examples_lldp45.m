% This script generates the Tables 3-5 of [1] illustrating the performance 
% of the codes LLDP1 and LLDP2 for different tolerances. 
%
% [1] Jacobian-free Locally Linearized Runge-Kutta method of 
%     Dormand and Prince for large systems of differential equations
%     by F.S. Naranjo-Noda and J.C. Jimenez
%

close all;
clear all;

abspath = which('run_examples_lldp45');
pos = strfind(abspath, filesep); pos = pos(end);
abspath = abspath(1:pos - 1);

cd(sprintf('%s%s%s',abspath,filesep,'..'));

s = [
 [abspath,filesep,'..',filesep,'demo_lld45;']...,
 [abspath,filesep,'..',filesep,'llint;']...,
 [abspath,filesep,'Brusselator;']...,
 [abspath,filesep,'Brusselator2D;']...,
 [abspath,filesep,'Burgers;']...,
 [abspath,filesep,'CUSP;']...,
 [abspath,filesep,'DND;']...,
 [abspath,filesep,'GrayScott2D;']...,
 [abspath,filesep,'utiles;']
];

path(s,path);
initllPaths(true);

cusp_example;
bg_example;
br2d_example;
gs2d_example;
br_example;
dnd_example;