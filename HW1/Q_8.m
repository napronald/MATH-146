clc 
clear all
close all

dt_01
dt_0125

function dt_01()
    disp('for dt = 0.1')
    dt = 0.1;
    t = 0; 
    Nsteps = 864000; 
    
    for j=1:Nsteps
        t = t + dt;
    end
    format long
    t;
    abs_err = abs(t-86400)
    rel_err = abs((t-86400)/86400)
end

% for dt = 0.125
function dt_0125()
    dt = 0.125; 
    t = 0;
    Nsteps = 86400/dt; 
    
    
    for j=1:Nsteps
        t = t + dt;
    end
    disp('for dt = 0.125')
    format long
    t;
    abs_err = abs(t-86400)
    rel_err = abs((t-86400)/86400)
end