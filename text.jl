t = Template(;
           user="amandanicotina",
           authors=["Amanda Nicotina", "David Dodel"],
           plugins=[
               License(name="MIT"),
               Git(),
               GitHubActions(),
           ],
       )



##########################################################################################
# Initial RF field
N   = 600;
αx  = π/2;
αy  = π/6;
t_c = 0.5;

time = range(0.0, t_c, N);
t    = time .- t_c/2;
rotx = rad2deg(αx)/360;
roty = rad2deg(αy)/360;
flip_x = rotx/diff(t)[1];
flip_y = roty/diff(t)[1];
BW_Hz = 300.0;
x     = BW_Hz.*t;
y     = BW_Hz.*t;
#B0    = [0.0];
B0    = [-150.0, -100.0, -50, 0.0, 50.0, 100.0, 150.0]; # [Hz]
Bz    = zeros(1,N);
B1x   = ((flip_x.*sinc.(x))./2π)'; # sinc(x) = sin(πx)/(πx)
B1y   = ((flip_y.*sinc.(y))./2π)';