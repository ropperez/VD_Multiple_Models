% This code provides the implementation for the VD-IMM and VD-GPB2, 
% developed in the paper "R. Pérez-Pérez, A. F. García-Fernández, "entitled 
% Variable Dimension IMM and GPB2 Filters for Tracking in Multiple Model" 
% accepted in IEEE Transactions on Aerospace and Electronic Systems. 
% Code author: R. Pérez-Pérez

clc
clear all      %#ok<CLALL>
close all

% Include the path
restoredefaultpath();
addpath(genpath('10-AUXILIARY'))
addpath(genpath('20-FILTER'))
addpath(genpath('30-GENERATION'))
addpath(genpath('40-MM'))
addpath(genpath('50-UTILS'))

% Simulation status
fprintf(['Simulation started at ' char(datetime('now')) '\n']);

% Enums
ENUM            = get_enum();
enum_filter     = ENUM.enum_filter;
enum_dyn_model  = ENUM.dynamic_model;

% Filters' tags
tag_filter      = fieldnames(enum_filter);

% Filter type simulation selection
sim_filter = [enum_filter.vd_imm enum_filter.vd_gpb2];

% Models to be generated
sim_r = [enum_dyn_model.cv enum_dyn_model.ct_r enum_dyn_model.ct_l];

% Run KF MM
VD_MM_simulation(sim_r,sim_filter,tag_filter,ENUM);

% Simulation status
fprintf(['Simulation ended at ' char(datetime('now')) '\n']);