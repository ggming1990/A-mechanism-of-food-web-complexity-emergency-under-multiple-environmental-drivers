function output=equation1(d,U1,U2,U3,ci,cj,ei,ej,H,u1,u2,FF) %d-food web matrix, U1-range of resource productivity, U2-range of ecosystem size, U3-range of disturbance, ci-basal species colonization rates, cj-consumers' colonization rates, ei-basal species mortality rate, ej-consumers' mortality rate, H-competition matrix, u1-predation rate of consumer preying on basal species, u2-predation rate between consumers, FF-save the outcome of interactive effects
C_length=length(U1);% length of resource productivity
l_con=length(U2);% length of ecosystem size
ll_length=length(U3);% length of disturbance
FDLmax=zeros(l_con,C_length); % maximum food chain length
d12=d; % food web matrix
[zwu,zwu1]=size(d);%total number of species
al=sum(d,1);
al1=sum(al==0);
M1=d(1:al1,al1+1:zwu);%interaction matrix between basal species and consumers
M2=d(al1+1:zwu,al1+1:zwu);%interaction matrix between consumers
[MM1,MM2]=size(M1);
m=MM1; % number of basal species
s1=MM2; % number of consumers
cc11 = repmat(ci,1,m); c111 = repmat(ci',m,1);
Tij=cc11.*H-c111.*H';
a2=rand(s1+m,1);%initial values
%%%%%%%%%%%%%%%%%%%%Interactive effects of S, R & D on FCL
for iii=1:l_con
    iii
    S=U2(iii); % ecosystem size
    for i111=1:C_length
        R=U1(i111); % resource productivity
        for i112=1:ll_length
            D=U3(i112); % disturbance
            rho0=(a2/sum(a2))*S; % initial species abundances
            t0 = 0:1000000 ; % total number of time steps
            [t,rho]=ode45(@(t,rho)equation(t,rho,S,Tij,ci,cj,ei,m,ej,M1,s1,M2,D,u1,u2,R),t0,rho0); %calling function "equation.m"
            qs=round(0.9*size(rho,1)); 
            qs21=rho(:,1:m);
            md21 = mean(qs21(qs:end,:)); % averaging the final 2000 time steps as the basal species abundances at steady state
            zw21=sum(md21>1e-6); % basal species with abundance less than 10-6 deemed as extinct
            md1 = mean(rho(qs:end,:)); % averaging the final 2000 time steps as the consumer species abundances at steady state
            zw1=find(md1>1e-6); % consumer species with abundance less than 10-6 deemed as extinct
           %Obtaining the interaction matrix between consumers
            AM = d12(zw1,:);
            d11 = AM(:,zw1);
            md12=md1(:,zw1);
            if  zw21==0 % if there is no basal species, then no species can survive
                FDLmaxSR(iii,i111)= 0;
                FDLmaxDS(i112,iii)= 0;
                FDLmaxRD(i111,i112)= 0;
            else
                md111=md12';%surviving species identity
                c=[1:zw21];%surviving basal sepcies
                d112=d11;%interaction matrix between consumers
                lp=size(d11,1);%number of surviving consumers
                trop=zeros(lp,1);%save trophic levels for species
                %setting basal species' trophic level as 1
                for is1=1:zw21
                    trop(is1,1)=1;
                end
                %estimating each species' trophic position
                for jj=1:10000
                    d112(:,c)=0;
                    ddd=d112(c,:);
                    ddd2=d112;
                    ddd2(c,:)=[];
                    if size(c,2)==1
                        ddds=ddd;
                    else
                        ddds=sum(ddd);
                    end
                    ddd2s=sum(ddd2);
                    dddpd=ddd2s./ddds;
                    cd=find(dddpd==0);
                    for imq=1:length(cd)
                        pm1= cd(imq);
                        pm2= find(d11(:,pm1)==1);
                        trop(pm1,1)=sum((md111(pm2,1)/sum(md111(pm2,1))).*trop(pm2,1))+1;
                    end
                    c=[c,cd];
                    CCC=size(c,2);
                    if CCC==size(d11,1)
                        break
                    end
                end
                zuida=find(sum(d11,2)==0);%Looking for top predators
                FDLmaxSR(iii,i111)= max(trop(zuida));%looking for the maximum trophic level under the interactive effect of ecosystem size and resource productivity
                FDLmaxDS(i112,iii)= max(trop(zuida));%looking for the maximum trophic level under the interactive effect of disturbance extent and ecosystem size
                FDLmaxRD(i111,i112)= max(trop(zuida));%looking for the maximum trophic level under the interactive effect of resource productivity and disturbance extent
            end
        end
    end
end
%%%%%%%%%save the outcome
if ll_length==1
    FF{1,1}=FDLmaxSR;
elseif l_con==1
    FF{2,1}=FDLmaxRD;
elseif C_length==1
    FF{3,1}=FDLmaxDS;
end
%%%%%%%%%%%%%%%%%%
output=FF;
