function Ldistr_5G_year1_solevacc() %was originally Year 0 but that's confusing

global Unvacc Left_C Doses Perc Left_budget Perc_totalpplnGr cov Bused

clear Doses
clear Perc
clear Left_budget  % left over budget, cash+doses translated in $
clear Left_C  % left over cash
clear Perc_totalpplnGr
clear Unvacc

scale=0.8;


%%Queue demand
%d1=

%dmd=[d1,d2,d3,d4,d5];

dmd=[0.2,0.2,0.2,0.2,0.2];                      % Case 1A - from NACI federal survey HZ coverage estimate
%dmd = [0.175,0.19825,0.2215,0.24475, 0.268]   % Case 1B - used equation of a line between (65, 0.175) and (69, 0.268) --> American data
%dmd =                                         % Case 1C - 
%dmd = 0.084 + (0.3-0.084)*rand(1,5);          % Case 1D - use 30% as a highball max coverage rate (US coverage goal for HZ)
overestB=0.5;
gamma1=1; %assume 100% monopoly of shingrix

%Gr is the matrix of [group sizes, demand in group, group label by age: 1= 60-61, 2=61-62 etc.]
Gr(1,:)=[152930,dmd(1),1]; % 65-66
Gr(2,:)=[148584,dmd(2),2]; %66-67
Gr(3,:)=[145840,dmd(3),3]; % 67-68
Gr(4,:)=[146319,dmd(4),4]; % etc
Gr(5,:)=[146037,dmd(5),5]; % etc
Gr;

l=1
for l=1:1

%We shuffle the order of the groups so we do not always start distributing to last group first.
j=1;
order=1:5;
order=order(randperm(length(order)));
for j=1:5
 % Gr(j,:)=Gr(order(j),:); % I am reindexing groups according to the shuffle of (1,2,3,4,5).
extract(j)=Gr(order(j),1)/1000;
demand(j)=Gr(order(j),2) ;
j=j+1;
end
extract;
demand;
%Gr


disp('Scaled group sizes after reorder!');
extract;

for i=1:5
    G(i)=extract(i)*demand(i);
    n(i)=floor(G(i));
    lambda(i)=floor((1-0.2*rand(1))*n(i));
    d(i)=(lambda(i)^(n(i))*exp(-lambda(i)))/factorial(n(i));
end
disp('Estimated demand in each group after reshuffle:');
G;
n
lambda
disp('queue demand in each group')
d

prob_outcome=[0.385;0.391;0.396;0.402;0.408]; %probab. of disease in group i leading to outcome "shingles" data from line btwn (50,0.3) and (85,0.5) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5724974/

%I compute the budget scaled in units of 1000, since i cut the groups by a factor of 1000.

% Overestimate budget at p1=300 for all;
p=scale*[300]; % price govt. pays = price of products on market * scale (scale = discount factor for gvt.)
av_price=gamma1*p(1);

  
 for i=1:5
     B(i)=overestB*(av_price*G(i));
 end
 
disp('Scaled budget in groups after reorder')
B
totB=sum(B);

%u units of resource = number of vaccine doses available as a Budget/price
%estimate in each group, of both kinds.

for i=1:5 
        u(i,1)=(gamma1*B(i))/p(1); % vaccines from producer 1 in group i
      % u(i,2)=(gamma2*B(i))/p(2); % vaccines from producer 2 in group i
end

% The vector w holds the scaled number of doses
disp('# of doses purchased with budget B for each group, of each vaccine');
u
for i=1:5
w_total(i,:)=min(floor(u(i,1)),lambda(i));
end
w_total

%%% Check budget computation:
for i=1:5
    chec(i)=dot(w_total(i,:),p');
end
chec
sum(chec);
B=chec

%Lines 76-80 tell me if I used less money than I allocated.
if (sum(chec) == totB)
  disp(' budget exhausted ')
  else
  disp('budget not exhausted')
  end
%%%%%%%

disp('efficacy of each vaccine in each age group');
e=[0.9;0.9;0.889;0.89;0.88];



%%------------------------------------------------------------
%% STAGES - we do the allocation in stages
%%------------------------------------------------------------
%
disp(' We initiate the number of doses we have to allocate');
W_initial=floor(u);

%% Starting at stage T = 5, since we're working with 5 groups. Note that we work with the reshuffled groups here.
for i=1:5
  disp('current group we allocate in')
 T=order(6-i)
 
%%revenue computation for current group T, in which we allocate:

disp(' Define the vector of resources compared with the size of current group:');
w=[(min(G(T),w_total(T,1)))];%,(min(G(T),w_total(T,2)))];

clear xt1 xt2

xt1=0:1:w(1);
%xt2=0:1:w(2);

xt=xt1; % this matrix holds the decision variables of the DP method; 

 clear rev  % rev = value of the objective function in current group T; we clear it for each group we allocate in, 
            %otherwise we bring over the updated revenue from the previous 
            % stage, which we don't want.

   for s=1:(w(1)+1)
   % for t=1:((w(2)+2)-s)
         rev(s,i)=(e(i,1)*xt1(s))*prob_outcome(i)*d(i);
    end 
  
   
disp('# of rows in rev= size of (xt1); # ofcols of rev = size of xt2')
sizext1=size(xt1);
%sizext2=size(xt2);
rev;

%si=size(rev(:,:,i));
revbis=rev;#reshape(rev(:,:,i),si(1),si(2));
 
 %%% Lines 132-136 compute the max values of the objective function in group T, and the allocations which achieve this value = 
 %i.e., we compute an optimal soln of the allocation problem
 disp('max value of ob. fction= M; optimal allocation values of xt1, xt2 in max_vacc')
 [M,I] = max(revbis(:));
 f(i)=M;
 [I_row, I_col] = ind2sub(size(revbis),I); 
 max_vacc1(T,1)=I_row-1;
% max_vacc2(T,1)=I_col-1;

 % resetting our resources: 

disp('Reset resources after every allocation in a group')
w_total(T,:)=max(W_initial(T,:)-[max_vacc1(T,1)],zeros(1,1));
clear revbis
 end
 
%% From here, we finish the assignments allocated to each group
 
disp('Final assignments per groups per vaccine')
 mv1=[max_vacc1]%,max_vacc2]
 
disp('Left over doses')
w_total


for k=1:5
    percc(k,:)=[floor(order(k)),w_total(k,:)];
    if (percc(k,1) == 1)
      Y(1,:)=percc(k,:);
      end
    if (percc(k,1) == 2)
      Y(2,:)=percc(k,:);
      end
      if (percc(k,1) == 3)
      Y(3,:)=percc(k,:);
      end
    if (percc(k,1) == 4)
      Y(4,:)=percc(k,:);
      end
      if (percc(k,1) == 5)
      Y(5,:)=percc(k,:);
      end
end

Doses(:,:,l)=Y(:,2:end);
clear percc
clear Y

disp('Percentage of people vaccinated in each group out who demanded it (the 20%)- reshuffled')
 for i=1:5
     % Percentage of people vaccinated in each group - assume likelihood to vacc if offered =100%:
      perc(i)= sum(mv1(i,:))/G(i);
      perc_totalpplnGr(i)= sum(mv1(i,:))/extract(i);
 end
 perc
 perc_totalpplnGr

for k=1:5
    percc(k,:)=[floor(order(k)),perc(k)];
    if (percc(k,1) == 1)
      Y(1,:)=percc(k,:);
      end
    if (percc(k,1) == 2)
      Y(2,:)=percc(k,:);
      end
      if (percc(k,1) == 3)
      Y(3,:)=percc(k,:);
      end
    if (percc(k,1) == 4)
      Y(4,:)=percc(k,:);
      end
      if (percc(k,1) == 5)
      Y(5,:)=percc(k,:);
      end
end


 Perc(:,l)=Y(:,2:end);
 clear percc
 clear Y

for k=1:5
    percc_tot(k,:)=[floor(order(k)),perc_totalpplnGr(k)];
    if (percc_tot(k,1) == 1)
      Y(1,:)=percc_tot(k,:);
      end
    if (percc_tot(k,1) == 2)
      Y(2,:)=percc_tot(k,:);
      end
      if (percc_tot(k,1) == 3)
      Y(3,:)=percc_tot(k,:);
      end
    if (percc_tot(k,1) == 4)
      Y(4,:)=percc_tot(k,:);
      end
      if (percc_tot(k,1) == 5)
      Y(5,:)=percc_tot(k,:);
      end
end
Perc_totalpplnGr(:,l)=Y(:,2:end);
 clear percc_tot
 clear Y
 
 disp('initial budget in groups')
 B
disp('Spent budget in reshuffled groups')
 for i=1:5
 Bopt(i)=mv1(i,1)*p(1);#+mv1(i,2)*p(2);
 Left(i)=B(i)-Bopt(i);
 Left_cash(i)=Left(i)-(w_total(i,1)*p(1));#+w_total(i,2)*p(2));
end

Bused(1:5,l)=Bopt';
 
clear w_total 
 for k=1:5
    percc1(k,:)=[floor(order(k)),Left(k)];
    if (percc1(k,1) == 1)
      Y1(1,:)=percc1(k,:);
      end
    if (percc1(k,1) == 2)
      Y1(2,:)=percc1(k,:);
      end
      if (percc1(k,1) == 3)
      Y1(3,:)=percc1(k,:);
      end
    if (percc1(k,1) == 4)
      Y1(4,:)=percc1(k,:);
      end
      if (percc1(k,1) == 5)
      Y1(5,:)=percc1(k,:);
      end
 end

Left_budget(:,l)=Y1(:,2:end);
 clear percc1 
clear Y1



 for k=1:5
    percc2(k,:)=[floor(order(k)),Left_cash(k)];
if (percc2(k,1) == 1)
      Y2(1,:)=percc2(k,:);
      end
    if (percc2(k,1) == 2)
      Y2(2,:)=percc2(k,:);
      end
      if (percc2(k,1) == 3)
      Y2(3,:)=percc2(k,:);
      end
    if (percc2(k,1) == 4)
      Y2(4,:)=percc2(k,:);
      end
      if (percc2(k,1) == 5)
      Y2(5,:)=percc2(k,:);
      end
  end

Left_C(:,l)=Y2(:,end);
clear percc2
clear Y2


end
 
 
%
%% -------------- End of STAGE Computations ------------
%
%%------------------------------------------------------------
%% PLOTS
%%------------------------------------------------------------
%
%% First, we reverse the shuffling (lines 184-206), so that we can plot the coverage for groups in proper age order.

Doses
Perc
Perc_totalpplnGr;
Left_budget;
Left_C;
%%%%

##########################G_young(1,:)=[157292,0.2,1]; % NEW AGE GROUP 64-65 FOR REDISTRIBUTION OF LEFTOVERS
##########################young_prob_outcome=0.38;
##########################e_young=[0.88,0.52];
##########################w=Doses(:,:,l)
##########################
##########################clear xt1 xt2
##########################
##########################xt1=0:1:w(1);
##########################xt2=0:1:w(2);
##########################
##########################xt=[xt1,xt2]; % this matrix holds the decision variables of the DP method; 
##########################
########################## for s=1:(w(1)+1)
##########################    for t=1:((w(2)+2)-s)
########################## rev(s,t,i)=(e_young(1)*xt1(s) + e_young(2)*xt2(t))*young_prob_outcome*0.2; %where 0.2 is demand (constant in this case)
##########################end 
########################## end
##########################  
##########################  si=size(rev(:,:,i));
##########################revbis=reshape(rev(:,:,i),si(1),si(2))
########################## 
########################## %%% Lines 132-136 compute the max values of the objective function in group T, and the allocations which achieve this value = 
########################## %i.e., we compute an optimal soln of the allocation problem
########################## disp('max value of ob. fction= M; optimal allocation values of xt1, xt2 in max_vacc')
########################## [M,I] = max(revbis(:));
########################## f(i)=M
########################## [I_row, I_col] = ind2sub(size(revbis),I); 
########################## max_young1=I_row-1
########################## max_young2=I_col-1
########################## cov(:,)=((max_young1+max_young2)*10^3)/G_young(1,1);
########################## % resetting our resources: 
##########################
##########################disp('Reset resources after every allocation in a group')
##########################Doses(:,:,l)=max(Doses(:,:,l)-[max_young1,max_young2],zeros(1,2))
##########################clear revbis



%%%%  
%%%%end
%%%%percc
%%%%Y
%%%%
%%%%%%%%%%% Data we need to roll over in the next year:%%%%%
#Unvacc_perc = 1- Perc;
#for k=1:5
#Unvacc(k,l)=Unvacc_perc(k,l)*Gr(k,1);
#end

disp('Unvaccinated numbers per (unshuffled) group after year 1')
for k=1:5
Unvacc(k,l)=(1-Perc_totalpplnGr(k,l))*Gr(k,1);
end
Unvacc;

%disp('Leftover budget from year 1')
%%Left
%%%
%%%%%%%Graphs in bar format, each group has its own bar- these are for percc.
%%%%%% 
v1=1:1:5;
 figure
  h=bar(1:1:5,Perc(:,l),"c");
  axis([0 6 0 1]);
 title('YEAR 1: Subgroup coverage as a % of total subgroup demand with 2 vaccines and random demand')
 figure
  h=bar(1:1:5,Perc_totalpplnGr(:,l),"r");
  axis([0 6 0 1]);
 title('YEAR 1: Subgroup coverage as a % of total subgroup population with 2 vaccines and random demand')
end  % end for the l  loop


function v=shuffle(v)
  v=v(randperm(length(v)))
end