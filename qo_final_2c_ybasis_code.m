%This is the code for the z basis so while I may write elements for
%the y components too, I will only end up using those for the z components
%we shall start by creating some variables
syms psi_f psi_i psi_n psi_dz psi_dy psi_d;
syms t dt t_rec step fns time;
syms dp p p_rec p_avg;
syms hbar gamma omega;
syms H_z H_y H_effz H_effy H;
syms dpl dm do dp_dag dm_dag do_dag P_m;
syms lp lm lo lp_dag lm_dag lo_dag;
syms U_z U_y;
syms eps1 eps2 eps3;
syms psq psq_avg msq var sdev sdev_tot;
syms size;

%now we shall create some of our matrices
dpl = (1/sqrt(2))*[0,0,0,0,0,0;
                  0,0,0,0,0,0;
                  0,0,0,0,0,0;
                  0,0,0,0,0,0;
                 -1,0,0,0,0,0;
                 0,-1,0,0,0,0];
dm = (1/sqrt(2))*[0,0,0,0,0,0;
                  0,0,0,0,0,0;
                  0,0,0,0,0,0;
                  0,1,0,0,0,0;
                  0,0,1,0,0,0;
                  0,0,0,0,0,0];
do = [0,0,0,0,0,0;
      0,0,0,0,0,0;
      0,0,0,0,0,0;
     -1,0,0,0,0,0;
      0,0,0,0,0,0;
      0,0,1,0,0,0];
dp_dag = dpl';
dm_dag = dm';
do_dag = do';
P_m = [1,0,0,0,0,0;
       0,1,0,0,0,0;
       0,0,1,0,0,0;
       0,0,0,0,0,0;
       0,0,0,0,0,0;
       0,0,0,0,0,0];

%now we shall input some of our constants
hbar = 1;
gamma = 1;
omega = 1;

%and now our hamiltonians 
H_z = -1*1i*hbar*omega*(1/sqrt(8))*(dpl - dp_dag + dm -dm_dag);
H_y = -1*hbar*omega*(1/2)*(do + do_dag);

%and from there our effective hamiltonains
H_effz = H_z -1*1i*hbar*gamma*(1/2)*P_m;
H_effy = H_y -1*1i*hbar*gamma*(1/2)*P_m;

%now we shall create our jump operators
lp = sqrt(gamma)*dpl;
lm = sqrt(gamma)*dm;
lo = sqrt(gamma)*do;
lp_dag = lp';
lm_dag = lm';
lo_dag = lo';

%adjusting our hamiltonian to the y basis

H = H_effy;

%now initializing our quantum sates
psi_i = (1/2)*[0;0;0;-1;sqrt(2);-1];
psi_f = (1/2)*[0;0;0;-1;sqrt(2);-1];
psi_n = (1/2)*[0;0;0;-1;sqrt(2);-1];

%adjusting step size
dt = .001;
t = 10;
step = t/dt;
fns = 100;

%then we create matrices to stor our probabilities
% and also a vector to represent a time scale
p_rec=zeros(step,fns);
t_rec=zeros(step,1);

%we create our time axis
for n=1:step
    t_rec(n) = (n-1)*dt;
end

%we initialize our probability
dp = 0;
p = 0;

%from there we get our time evolutino matrix approximation
U = eye(6)-1i*(1/hbar)*dt*H;

%next we initialize our random numbers
eps1 = 0;
eps2 = 0;


%next we create our dark state eigenvectors to compare against;
psi_dz = (1/sqrt(2))*[0,0,0,1,0,1];
psi_dy = [0,0,0,0,1,0];

%choosing the correct coordinate system
psi_d = psi_dy;


%now we must create multiple wave functions
for m = 1:fns
   %within this we must create a single quantum trajectory
   %and we must reset for each one
   psi_i = (1/2)*[0;0;0;-1;sqrt(2);-1];
   psi_f = (1/2)*[0;0;0;-1;sqrt(2);-1];
   psi_n = (1/2)*[0;0;0;-1;sqrt(2);-1];
   dp = 0;
   p = 0;
   eps1 = rand;
   time = 0.0001;
   size = 0;
   eps2 = 0;
   eps3 = 0;

   %create an index to allow for p to be inserted into its tracking matrix
   n=1;
   
   %create a while loop so that this runs during our time interval
   
   while time < t-dt
   
       %set the next wavefunction to be the same as the one at the end of
       %the previous steo
       psi_i = psi_f;
       %measure the size of the wavefunction
       size = psi_i'*(U'*U)*psi_i;
       %create the if statement as to whether a jump occurs or not
       if size > eps1
           psi_f = U*psi_i;
       elseif size < eps1
          eps2 = rand;
          %this part we need to equally weight and apply each jump operator, this relies on the chosen coordinate system
          eps3 = 4*eps2;
          %now we create an equal liklihood of every jump operator
          if eps3 < 1
            psi_f = lp*psi_i;
            if norm(psi_f) == 0
                psi_f = U*psi_i;
            else 
                psi_f = lp*psi_i;
                psi_f = psi_f/norm(psi_f);
            end
            eps1 = rand;
              
              
          elseif eps3 > 1 && eps3 < 2
            psi_f = lo*psi_i;
            if norm(psi_f) == 0
                psi_f = U*psi_i;
            else 
                psi_f = lo*psi_i;
                psi_f = psi_f/norm(psi_f);
            end
            eps1 = rand;
              
              
          elseif eps3 > 2 && eps3 < 3
            psi_f = lm*psi_i;
            if norm(psi_f) == 0
                psi_f = U*psi_i;
            else 
                psi_f = lm*psi_i;
                psi_f = psi_f/norm(psi_f);
            end
            eps1 = rand;
              
              
          elseif eps3 > 3 
            psi_f = lo_dag*psi_i;
             if norm(psi_f) == 0
                psi_f = U*psi_i;
            else 
                psi_f = lo_dag*psi_i;
                psi_f = psi_f/norm(psi_f);
            end
            eps1 = rand;
              
              
          
          end
       end
       
   %now we measure the probability of the excited state
   psi_n= psi_f/norm(psi_f);
   p = (abs(psi_d*psi_n))^(2);
   p_rec(n,m) = p;
   %and advance time and the tracker
   n=n+1;
   time = time + dt; 
   
   end
end

p_avg=(1/fns)*sum(p_rec,2);

plot(t_rec,p_avg);

title('reproduction of Molmer Fig.11 one hundred trajectory y basis');
xlabel('time(unit:1/gamma)')
ylabel('Probability of the dark state');
