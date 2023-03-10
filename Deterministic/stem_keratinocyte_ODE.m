kp = 0.0385;
kp_dis = 0.0165;
kd = 0.0165;
kd_dis = 0.0385;
ka = 0.195;
kc = 0.00055;
ks = 0.065;
initial = [40,2700];
initial_dis = [40, 2679]
tspan = [0 260];
tdis = [50 60];

[t,x] = ode45(@(t,x) sk(kp, kd, kc, ka, ks, x, t), tspan, initial);
[t_dis, x_dis] = ode45(@(t,x) sk(kp_dis, kd_dis, kc, ka, ks, x, t), tdis, initial_dis)

plot(t_dis, x_dis(:,1))

yyaxis left
plot(t, x(:,1))
ylim([0 60])
ylabel("Stem Cell Population")
yyaxis right
plot(t, x(:,2))
ylim([0 3500])
ylabel("Keratinocyte Population")
%% functions
function dxdt = sk(kp, kd, kc, ka, ks, x, t);
    ds_dt = kp*x(1) - kd*x(1) - kc * x(1)^2;
    dk_dt = 32*kd*x(1) + 16*ka*x(1) + 32*kc*(x(1)^2) - ks*x(2);

    dxdt = [ds_dt;dk_dt];
end