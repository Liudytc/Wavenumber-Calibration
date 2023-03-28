clear;close all; 
load rawdata.mat; %rawdata=100 spectrum 
peak_pixel_approx_position=[0;19;42;49;75;100;111;131;156;168;181;197;206;232;243;260;532;561;572;625]';
peak_numbers=size(peak_pixel_approx_position,2);
for i=1:100
    temp_point=fpeak(1:1024,rawdata(i,:),3,[1,600,6000,10e4]);
    first_point=temp_point(1,1);
    peak_range1(:,:)=[first_point+peak_pixel_approx_position-3;first_point+peak_pixel_approx_position-2;first_point+peak_pixel_approx_position-1;first_point+peak_pixel_approx_position;first_point+peak_pixel_approx_position+1;first_point+peak_pixel_approx_position+2;first_point+peak_pixel_approx_position+3];
   peak_range=peak_range1';
    for j=1:peak_numbers
    peak_pixel_position{i}(j,:)=fpeak(1:1024,rawdata(i,:),3,[peak_range(j,1),peak_range(j,7),1,10e4]);    
    end
    peak_positions_pixel(i,:)=peak_pixel_position{i}(:,1); %Find the real peak pixels
end
%% Find the accurate peak position
% chosen_k=K_chosen(peak_positions_pixel,peak_numbers,rawdata);%Find the suitable parameters for lorentzfit
load ace_g300_k_v1.mat;
wind1=[3,3,3,3,3,3,3,2,2,3,3,3,3,3,3,3,3,3,2,3]; % The left part of window
wind2=[3,3,3,3,3,3,3,2,2,3,3,3,3,3,3,3,3,3,3,3]; % The right part of window
for i=1:100
    for j=1:peak_numbers
    p=peak_positions_pixel(i,j)-wind1(j);
    q=peak_positions_pixel(i,j)+wind2(j);    
    y=rawdata(i,p:q);   
    x=p:1:q;
    P2=peak_positions_pixel(i,j);P3=chosen_k(i,j);C=0;
    [height,pos]= max(y);
    P1=P3*height;
    P0=[P1,P2,P3,C];
    [yprime2 P resnorm2 residual2] = lorentzfit(x,y,P0);
    [peak_positions_subpixel(i,j)]=P(2);  
    end
end
peak_positions_subpixel=peak_positions_subpixel-512; 
%% Parameters 
wavenumber_NIST=([213.3,329.2,465.1,504.0,651.6,797.2,857.9,968.7,1105.5,1168.5,1236.8,1323.9,1371.5,1515.1,1561.5,1648.4,2931.1,3064.6,3102.4,3326.6]);
laser=532e-9;
current_grating=300;
d=1e-3/current_grating;
T = 26.00e-6;
fy_centre = 21.88;
alfa=fy_centre/2;
rotate_centre=5;
rotate_range = 4;
rotate_n=501;
rotate_step_size = rotate_range/(rotate_n-1);
c_centre = 0;
c_range = 200;
c_n = 1001;
c_step_size = c_range/(c_n-1);
f=0.5;
%% First guess
l=peak_numbers/2; %Choose the middle peak
rotate=rotate_centre+rotate_range/2;
tenth_peak_NIST =-(f/T)*tand(-alfa + rotate+asind(-laser./(d-1e2.*wavenumber_NIST(l).*laser.*d)-sind(-rotate-alfa)))+c_centre;
for i=1:100
MIN(i) = sum(abs(peak_positions_subpixel(i,l)-tenth_peak_NIST));
    for rotate=rotate_centre-rotate_range/2:rotate_step_size:rotate_centre+rotate_range/2
        for c=c_centre-c_range/2:c_step_size:c_centre+c_range/2            
        tenth_peak_NIST =-(f/T)*tand(-alfa + rotate+asind(-laser./(d-1e2.*wavenumber_NIST(l).*laser.*d)-sind(-rotate-alfa)))+c; 
            if sum(abs(peak_positions_subpixel(i,l)-tenth_peak_NIST))<MIN(i)
                MIN(i) = sum(abs(peak_positions_subpixel(i,l)-tenth_peak_NIST));
                temp=tenth_peak_NIST;
                C(i) = c; 
                ROTATE(i)=rotate;
            end
        end
    end  
end
%% REFINED SEARCH 
laser_centre = 532e-9;
laser_range =0e-9;
laser_n = 1;
laser_step_size = 1;

rotate_centre = round(ROTATE.*1e4)./1e4; %TAKEN FROM FIRST GUESS
rotate_range =0.001;
rotate_n=101;
rotate_step_size =rotate_range/(rotate_n-1);
c=C;
d = 1e-3./current_grating;
N=1;M=1;L=peak_numbers;c=C;
for i=1:100
MIN1{i}=ones(L,1)*1e3;MIN2{i}=ones(L,1)*1e3;MIN3{i}=ones(L,1)*1e3;
peak_positions_subpixel_NIST_temp{i}=ones(L,L-1);
First{i}=ones(L,2);Second{i}=ones(L,3);Third{i}=ones(L,4);
First_inv{i}=ones(L,2);Second_inv{i}=ones(L,3);Third_inv{i}=ones(L,4);
first_pixel{i}=ones(L,L-1);second_pixel{i}=ones(L,L-1);third_pixel{i}=ones(L,L-1);
end
ROTATE1=zeros(100,L);ALFA1=zeros(100,L); D1=zeros(100,L);  
ROTATE2=zeros(100,L);ALFA2=zeros(100,L); D2=zeros(100,L); 
ROTATE3=zeros(100,L);ALFA3=zeros(100,L); D3=zeros(100,L);
fun1=@(K)wavenumber_NIST(K);   
fun2=@(i,K)peak_positions_subpixel(i,K);

parfor i=1:100
        for P=N:L
        K=(N:M:L);    
        K(P)=[];  
        for rotate=rotate_centre(i)-rotate_range/2:rotate_step_size:rotate_centre(i)+rotate_range/2
            peak_positions_subpixel_NIST_temp{i}(P,:) =-(f/T)*tand(-alfa + rotate+asind(-laser./(d-1e2.*fun1(K).*laser.*d)-sind(-rotate-alfa)))+c(i);
            first=polyfit(peak_positions_subpixel_NIST_temp{i}(P,:),fun2(i,K),1);
            first_inv=[1/first(1),-first(2)/first(1)];  
            X_first=polyval(first_inv,fun2(i,K));
            temp_wavenum1 = (1./(d./laser.*(sind(atand((c(i).*T-X_first.*T)./f)+alfa-rotate)+sind(-alfa-rotate)))+1).*1e-2./laser;
            Er1_inv=mean(abs(temp_wavenum1-fun1(K)));
            Er1_raw=abs(temp_wavenum1-fun1(K));
            first_all=polyval(first,peak_positions_subpixel_NIST_temp{i}(P,:));            
            if Er1_inv<MIN1{i}(P,:)
                 MIN1{i}(P,:) =Er1_inv;
                 ROTATE1(i,P)=rotate;%
                 ALFA1(i,P)=alfa;%
                 D1(i,P)=d;
                 Laser1(i,P)=laser;
                 First{i}(P,:)=first;
                 First_inv{i}(P,:)=first_inv;
                 first_pixel{i}(P,:)=first_all;
                 Error1_raw{i}(P,:)=Er1_raw;
            end
        end
        end
end

for i=1:100
        for P=N:L    
           first_subpixel_polyed(i,P)=polyval(First_inv{i}(P,:),peak_positions_subpixel(i,P));
           our_calibrated_wavenum1(i,P) = (1./(D1(i,P)./Laser1(i,P).*(sind(atand((c(i).*T-first_subpixel_polyed(i,P).*T)./f)+ALFA1(i,P)-ROTATE1(i,P))+sind(-ALFA1(i,P)-ROTATE1(i,P))))+1).*1e-2./Laser1(i,P);
           error1_original1(i,P) = our_calibrated_wavenum1(i,P) - wavenumber_NIST(P);
        end
end
%% traditional
for i=1:100
    for P=N:L
    K=(N:M:L); K(P)=[]; 
        temp1=polyfit(peak_positions_subpixel(i,K),wavenumber_NIST(K),1);
        temp2(i,P)=polyval(temp1,peak_positions_subpixel(i,P));
        error2(i,P) =mean(abs(temp2(i,P)-wavenumber_NIST(P))); 

        temp3=polyfit(peak_positions_subpixel(i,K),wavenumber_NIST(K),2);
        temp4(i,P)=polyval(temp3,peak_positions_subpixel(i,P));    
        error3(i,P) =mean(abs(temp4(i,P)-wavenumber_NIST(P))); 

        temp5=polyfit(peak_positions_subpixel(i,K),wavenumber_NIST(K),3);
        temp6(i,P)=polyval(temp5,peak_positions_subpixel(i,P));
        error4(i,P) =mean(abs(temp6(i,P)-wavenumber_NIST(P))); 

        temp7=polyfit(peak_positions_subpixel(i,K),wavenumber_NIST(K),4);
        temp8(i,P)=polyval(temp7,peak_positions_subpixel(i,P));
        error5(i,P) =mean(abs(temp8(i,P)-wavenumber_NIST(P))); 
        
        temp9=polyfit(peak_positions_subpixel(i,K),wavenumber_NIST(K),5);
        temp10(i,P)=polyval(temp9,peak_positions_subpixel(i,P));
        error6(i,P) =mean(abs(temp10(i,P)-wavenumber_NIST(P))); 
        
        temp11=polyfit(peak_positions_subpixel(i,K),wavenumber_NIST(K),6);
        temp12(i,P)=polyval(temp11,peak_positions_subpixel(i,P));
        error7(i,P) =mean(abs(temp12(i,P)-wavenumber_NIST(P)));
        
        temp13=polyfit(peak_positions_subpixel(i,K),wavenumber_NIST(K),7);
        temp14(i,P)=polyval(temp13,peak_positions_subpixel(i,P));
        error8(i,P) =mean(abs(temp14(i,P)-wavenumber_NIST(P)));
    end
end
first_result=mean(mean(abs(error1_original1)));
order1=mean(mean(abs(error2)));order2=mean(mean(abs(error3)));order3=mean(mean(abs(error4)));
order4=mean(mean(abs(error5)));order5=mean(mean(abs(error6)));order6=mean(mean(abs(error7)));order7=mean(mean(abs(error8)));
for i=1:100
std1(i)=std(error1_original1(i,:));    
std_orignial1(i)=std(temp2(i,:)-wavenumber_NIST);
std_orignial2(i)=std(temp4(i,:)-wavenumber_NIST);
std_orignial3(i)=std(temp6(i,:)-wavenumber_NIST);
std_orignial4(i)=std(temp8(i,:)-wavenumber_NIST);
std_orignial5(i)=std(temp10(i,:)-wavenumber_NIST);
std_orignial6(i)=std(temp12(i,:)-wavenumber_NIST);
std_orignial7(i)=std(temp14(i,:)-wavenumber_NIST);
end
STD_e1=mean(std1);STD_o1=mean(std_orignial1);STD_o2=mean(std_orignial2);STD_o3=mean(std_orignial3);
STD_o4=mean(std_orignial4);STD_o5=mean(std_orignial5);STD_o6=mean(std_orignial6);STD_o7=mean(std_orignial7);
for i=1:100
rmse1(i)=sqrt(mean(error1_original1(i,:).^2));    
rmse_orignial1(i)=sqrt(mean((temp2(i,:)-wavenumber_NIST).^2));
rmse_orignial2(i)=sqrt(mean((temp4(i,:)-wavenumber_NIST).^2));
rmse_orignial3(i)=sqrt(mean((temp6(i,:)-wavenumber_NIST).^2));
rmse_orignial4(i)=sqrt(mean((temp8(i,:)-wavenumber_NIST).^2));
rmse_orignial5(i)=sqrt(mean((temp10(i,:)-wavenumber_NIST).^2));
rmse_orignial6(i)=sqrt(mean((temp12(i,:)-wavenumber_NIST).^2));
rmse_orignial7(i)=sqrt(mean((temp14(i,:)-wavenumber_NIST).^2));
end
RMSE_e1=mean(rmse1);RMSE_o1=mean(rmse_orignial1);RMSE_o2=mean(rmse_orignial2);RMSE_o3=mean(rmse_orignial3);
RMSE_o4=mean(rmse_orignial4);RMSE_o5=mean(rmse_orignial5);RMSE_o6=mean(rmse_orignial6);RMSE_o7=mean(rmse_orignial7);
%% result=[Our method, 1st order, 2nd order,3rd order,4th order,5th order,6th order,7th order]
RMSE_ace_g300_loo_peaks=[RMSE_e1,RMSE_o1,RMSE_o2,RMSE_o3,RMSE_o4,RMSE_o5,RMSE_o6,RMSE_o7];
MAE_ace_g300_loo_peaks=[first_result,order1,order2,order3,order4,order5,order6,order7];
STD_ace_g300_loo_peaks=[STD_e1,STD_o1,STD_o2,STD_o3,STD_o4,STD_o5,STD_o6,STD_o7];

