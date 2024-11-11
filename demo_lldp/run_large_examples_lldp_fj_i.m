% This script generates the Tables for illustrating the performance 
% of the codes LLDP1 and LLDP2 for different tolerances. 
%
% [1] New class of Jacobian-free high order local
% Linearization methods for differential equations
% by F.S. Naranjo-Noda and J.C. Jimenez
%

close all;
clear all;

abspath = which('run_large_examples_lldp_fj_i');
pos = strfind(abspath, filesep); pos = pos(end);
abspath = abspath(1:pos - 1);

cd(sprintf('%s%s%s',abspath,filesep,'..'));

s = [
 [abspath,filesep,'..',filesep,'demo_lldp;']...,
 [abspath,filesep,'..',filesep,'llint;']...,
 [abspath,filesep,'Brusselator2D;']...,
 [abspath,filesep,'GrayScott2D;']...,
 [abspath,filesep,'utiles;']
];

path(s,path);
initllPaths(true);

br2d_example_large;
gs2d_example_large;
