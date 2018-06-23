M=100;                 % number of synapses per neuron
D=5;                   % maximal conduction delay 
% excitatory neurons   % inhibitory neurons      % total number 
Ne=800;                Ni=200;                   N=Ne+Ni;
a=[0.02*ones(Ne,1);    0.1*ones(Ni,1)];
d=[   8*ones(Ne,1);    2*ones(Ni,1)];
sm=4;                 % maximal synaptic strength

post=ceil([N*rand(Ne,M);Ne*rand(Ni,M)]); 
%s=[ones(Ne,M);-ones(Ni,M)];         % synaptic weights
s=[1*ones(Ne,M);-1*ones(Ni,M)];
sd=zeros(N,M);                      % their derivatives

%{
for i=1:N
  if i<=Ne
    for j=1:D
      delays{i,j}=M/D*(j-1)+(1:M/D);
    end;
  else
    delays{i,1}=1:M;
  end;
  pre{i}=find(post==i&s>0);             % pre excitatory neurons
  aux{i}=N*(D-1-ceil(ceil(pre{i}/N)/(M/D)))+1+mod(pre{i}-1,N);
end;
%}

delays = cell(N,D);
for i=1:Ne
    p=randperm(N);
    post(i,:)=p(1:M);
    for j=1:M
        delays{i, ceil(D*rand)}(end+1) = j;  % Assign random exc delays
    end;
end;
for i=Ne+1:N
    p=randperm(Ne);
    post(i,:)=p(1:M);
    delays{i,1}=1:M;                    % all inh delays are 1 ms.
end;

pre = cell(N,1);
aux = cell(N,1);
for i=1:Ne
    for j=1:D
        for k=1:length(delays{i,j})
            pre{post(i, delays{i, j}(k))}(end+1) = N*(delays{i, j}(k)-1)+i;
            aux{post(i, delays{i, j}(k))}(end+1) = N*(D-1-j)+i; % takes into account delay
        end
    end
end
  
STDP = zeros(N,1001+D);
v = -65*ones(N,1);                      % initial values
u = 0.2.*v;                             % initial values
firings=[-D 0];                         % spike timings

%---------------
% new stuff related to DA-STDP
T=3600*10;         % the duration of experiment
DA=0;           % level of dopamine above the baseline

rew=[];
neg_rew=[];

%selecting excitatory neurons for action representation
action_neurons= randperm(Ne,300);
n1_right=action_neurons(1:50);
n1_left=action_neurons(51:100);
%n1_stim_1=action_neurons(101:150);
%n1_stim_2=action_neurons(151:200);
%n1_stim_3=action_neurons(201:250);
%n1_stim_4=action_neurons(251:300);

count_left=0;
count_right=0;
count_no_action=0;
incorrect_L=0;
incorrect_R=0;

dlmwrite('action_neurons.csv',action_neurons,'delimiter',',','-append');

%4 inputs representation
I_input=14*(ones(N,4));
for i=1:1000
  if ~any(i==action_neurons(101:150))
    I_input(i,1)=0; %earlier kept 0
  end
  if ~any(i==action_neurons(151:200))
    I_input(i,2)=0;
  end
   if ~any(i==action_neurons(201:250))
    I_input(i,3)=0;
   end
   if ~any(i==action_neurons(251:300))
    I_input(i,4)=0;
   end
end

%shist_R=zeros(1000*T,2);
%shist_L=zeros(1000*T,2);
stimulus_time=0;


%--------------
alt=0;
input=1;

window=40;

for sec=1:T  % simulation of 1 hour 
     fired_R=0;
     fired_L=0;
     avg_firing=0; 
  for t=1:1000                          % simulation of 1 sec
        if(mod(sec,5)==0)
            if(t==1)
                if(alt==0)
                    input=1;
                    alt=1;
                else
                    input=2;
                    alt=0;
                end
                input2= randperm(2,1)+2;
                I=max(I_input(:,input),I_input(:,input2));
                I=I_input(:,input);
                stimulus_time=sec*1000+t;
                fprintf('sec:%d\t correct L:%d \tR:%d \tNA:%d \n%d\t Incorrect L:%d R:%d \n',sec,count_left,count_right,count_no_action,input,incorrect_L,incorrect_R);
            else
                if(t>=window+50)
                   I=15*(rand(N,1)-0.5);
                end
            end
        else
            I=15*(rand(N,1)-0.5);
        end
        
    fired = find(v>=30);                % indices of fired neuron
    v(fired)=-65;  
    u(fired)=u(fired)+d(fired);
    STDP(fired,t+D)=0.1;
    for k=1:length(fired)
      sd(pre{fired(k)})=sd(pre{fired(k)})+STDP(N*t+aux{fired(k)});
    end
    firings=[firings;t*ones(length(fired),1),fired];
    k=size(firings,1);
    while firings(k,1)>t-D
      del=delays{firings(k,2),t-firings(k,1)+1};
      ind = post(firings(k,2),del);
      I(ind)=I(ind)+s(firings(k,2), del)';
      sd(firings(k,2),del)=sd(firings(k,2),del)-1.2*STDP(ind,t+D)';
      k=k-1;
    end
    
    v=v+0.5*((0.04*v+5).*v+140-u+ I);    % for numerical 
    %v=v+0.4*((0.04*v+5).*v+140-u+ I);    % stability time 
    u=u+a.*(0.2*v-u);                   % step is 0.5 ms
    STDP(:,t+D+1)=0.95*STDP(:,t+D);     % tau = 20 ms
    
    if(t<=window+20 && mod(sec,5)==0)
        temp_R=0;
        temp_L=0;
        avg_firing=avg_firing+size(fired,1);
        for i=1:size(fired,1)
            if any(fired(i)==n1_right)
                temp_R=temp_R+1;
            end
            if any(fired(i)==n1_left)
                temp_L=temp_L+1;
            end
        end
        
        fired_R=fired_R+temp_R;
        fired_L=fired_L+temp_L;
        
        if (t==window+20)  %maybe take a extended window for counting the no of spikes in action neurons 
            action_executed=[0,0,0]; %[left,right,no action]
            correct=0; %correct=1 for correct case and correct=0 for wrong case  
            avg_firing=avg_firing/50;
            fprintf('size: avg: %d fired_ L: %d, R: %d \n',avg_firing,fired_L,fired_R);
            reward_time=-999;
                             
            if(fired_L>fired_R)
                action_executed=[1,0,0];
                reward_time=sec*1000+t+floor(fired_R/fired_L*1000);
                if (input==1) 
                    count_left=count_left+1;
                    correct=1;
                    %fprintf('action:left\tcount:%d \n',count_left);
                    %reward
                    
                    rew=[rew,reward_time];   %delay of 1 second
                else
                    incorrect_L=incorrect_L+1;
                    neg_rew=[neg_rew,reward_time];
                end
    
            elseif(fired_R>fired_L)
                action_executed=[0,1,0];
                reward_time=sec*1000+t+floor(fired_L/fired_R*1000);
                if (input==2)
                  correct=1;
                  count_right=count_right+1;
                    %fprintf('action:right \t  count:%d \n',count_right);
                    %reward
                  
                    rew=[rew,reward_time];   %delay of 1 second
                else
                     incorrect_R=incorrect_R+1;
                     neg_rew=[neg_rew,reward_time];
                end
            else 
                 action_executed=[0,0,1];
                 count_no_action=count_no_action+1;
                %fprintf('action:NA \t count:%d \n',count_no_action);
            end
            record=[stimulus_time,reward_time,correct,input,action_executed,avg_firing,fired_L,fired_R,input2];
            dlmwrite('data.csv',record,'delimiter',',','-append');
        end
    end
  
  if any(rew==sec*1000+t)
      DA=DA+0.5;
  end
  
  if any(neg_rew==sec*1000+t)
    DA=DA-0.5; 
  end
  
%{    
  if any(neg_rew==sec*1000+t)
      DA=DA*0.95;
    weight_record=[sec*1000+t,size(find(s>0)),s_sum_L,s_sum_R,sd_sum_L,sd_sum_R,DA];
    dlmwrite('weight.csv',weight_record,'delimiter',',','-append'); 
  end
%}  
  if (mod(t,10)==0)
       % s(1:Ne,:)=max(0,min(sm,s(1:Ne,:)+(0.002+DA)*sd(1:Ne,:)));
        s(1:Ne,:)=max(0,s(1:Ne,:)+(0.002+DA)*sd(1:Ne,:));
        sd=0.9*sd;
  end
  DA=DA*0.993;
%{  
    s_sum_R=0;
    s_sum_L=0;
    sd_sum_R=0;
    sd_sum_L=0;
 
    for k=1:50
        for j=1:100
            s_sum_L=s_sum_L+s(n1_left(k),j);
            sd_sum_L=sd_sum_L+sd(n1_left(k),j);
            
            s_sum_R=s_sum_R+s(n1_right(k),j);
            sd_sum_R=sd_sum_R+sd(n1_right(k),j);
        end
    end
    
    s_sum_L=s_sum_L/(100*50);
    s_sum_R=s_sum_R/(100*50);
    sd_sum_R=sd_sum_R/(100*50);
    sd_sum_L=sd_sum_L/(100*50);
    
    weight_record=[sec*1000+t,size(find(s>0)),s_sum_L,s_sum_R,sd_sum_L,sd_sum_R,DA];  
    dlmwrite('weight.csv',weight_record,'delimiter',',','-append');
    
    shist_R(sec*1000+t,:)=[s_sum_R,sd_sum_R];
    shist_L(sec*1000+t,:)=[s_sum_L,sd_sum_L];
%}
  end
%{
% ---- plot -------
  subplot(2,2,1)
  plot(firings(:,1),firings(:,2),'.');
  axis([0 1000 0 N]); 
  subplot(2,2,2)
  plot(0.001*(1:(sec*1000+t)),DA, 0.001*rew,0*rew,'rx');
  subplot(2,2,3);
  plot(0.001*(1:(sec*1000+t)),shist_R(1:sec*1000+t,:),0.001*(1:(sec*1000+t)),shist_L(1:sec*1000+t,:), 0.001*rew,0*rew,'rx',0.001*neg_rew,0*neg_rew,'bx');
  subplot(2,2,4);
  hist(s(find(s>0)),sm*(0.01:0.01:1)); % only excitatory synapses
  drawnow;
%}
% ---- end plot ------

  STDP(:,1:D+1)=STDP(:,1001:1001+D);
  ind = find(firings(:,1) > 1001-D);
  firings=[-D 0;firings(ind,1)-1000,firings(ind,2)];
end
