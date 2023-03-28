clear;close all; 
load rawdata.mat; %rawdata=100 spectrum 
test =[11:20]; %test
train=[1:10]; %training
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
d = 1e-3./current_grating;
laser_centre = laser;
laser_range =0e-9;
laser_n = 1;
laser_step_size = 1;

rotate_centre = round(ROTATE.*1e4)./1e4; %TAKEN FROM FIRST GUESS
rotate_range =0.001;
rotate_n=201;
rotate_step_size =rotate_range/(rotate_n-1);
c=C;%TAKEN FROM FIRST GUESS

MIN1(1:100)=1e3;MIN2(1:100)=1e3;MIN3(1:100)=1e3;
peak_positions_subpixel_NIST_temp=zeros(100,size(test,2));
ROTATE1=zeros(1,100);ALFA1=zeros(1,100); D1=zeros(1,100); 
ROTATE2=zeros(1,100);ALFA2=zeros(1,100); D2=zeros(1,100); 
ROTATE3=zeros(1,100);ALFA3=zeros(1,100); D3=zeros(1,100); 
First=zeros(100,2);first_pixel=zeros(100,size(test,2));First_inv=zeros(100,2);
Second=zeros(100,3);second_pixel=zeros(100,size(test,2));Second_inv=zeros(100,3);
Third=zeros(100,4);third_pixel=zeros(100,size(test,2));Third_inv=zeros(100,4);
%%
parfor i=1:100
    for laser=laser_centre-laser_range/2:laser_step_size:laser_centre+laser_range/2
        for rotate=rotate_centre(i)-rotate_range/2:rotate_step_size:rotate_centre(i)+rotate_range/2
        peak_positions_subpixel_NIST_temp(i,:) =-(f/T)*tand(-alfa + rotate+asind(-laser./(d-1e2.*wavenumber_NIST(test).*laser.*d)-sind(-rotate-alfa)))+c(i);
        first=polyfit(peak_positions_subpixel_NIST_temp(i,:),peak_positions_subpixel(i,test),1);
        first_inv=[1/first(1),-first(2)/first(1)];
        X_first=polyval(first_inv,peak_positions_subpixel(i,test));
        temp_wavenum1 = (1./(d./laser.*(sind(atand((c(i).*T-X_first.*T)./f)+alfa-rotate)+sind(-alfa-rotate)))+1).*1e-2./laser;
        Er1_inv=mean(abs(temp_wavenum1-wavenumber_NIST(test)));
        first_all=polyval(first,peak_positions_subpixel_NIST_temp(i,:));    
        if Er1_inv<MIN1(i)
        MIN1(i) =Er1_inv;
        ROTATE1(i)=rotate;%
        ALFA1(i)=alfa;%
        D1(i)=d;
        Laser1(i)=laser;
        First(i,:)=first;
        First_inv(i,:)=first_inv;
        first_pixel(i,:)=first_all;
        end
        end
    end
end
for i=1:100
X1=polyval(First_inv(i,:),peak_positions_subpixel(i,train));
wavenumber1(i,:)=(1./(D1(i)./Laser1(i).*(sind(atand((c(i).*T-X1.*T)./f)+ALFA1(i)-ROTATE1(i))+sind(-ALFA1(i)-ROTATE1(i))))+1).*1e-2./Laser1(i);
error1_original1(i,:)=mean(sqrt((wavenumber1(i,:)-wavenumber_NIST(train)).^2));
end

%% traditional
for i=1:100
        temp1=polyfit(peak_positions_subpixel(i,test),wavenumber_NIST(test),1);
        temp2(i,:)=polyval(temp1,peak_positions_subpixel(i,train));
        error2(i,:) =mean(abs(temp2(i,:)-wavenumber_NIST(train))); 

        temp3=polyfit(peak_positions_subpixel(i,test),wavenumber_NIST(test),2);
        temp4(i,:)=polyval(temp3,peak_positions_subpixel(i,train));    
        error3(i,:) =mean(abs(temp4(i,:)-wavenumber_NIST(train)));

        temp5=polyfit(peak_positions_subpixel(i,test),wavenumber_NIST(test),3);
        temp6(i,:)=polyval(temp5,peak_positions_subpixel(i,train));
        error4(i,:) =mean(abs(temp6(i,:)-wavenumber_NIST(train)));
        
        temp7=polyfit(peak_positions_subpixel(i,test),wavenumber_NIST(test),4);
        temp8(i,:)=polyval(temp7,peak_positions_subpixel(i,train));
        error5(i,:) =mean(abs((temp8(i,:)-wavenumber_NIST(train))));
        
        temp9=polyfit(peak_positions_subpixel(i,test),wavenumber_NIST(test),5);
        temp10(i,:)=polyval(temp9,peak_positions_subpixel(i,train));
        error6(i,:) =mean(abs(temp10(i,:)-wavenumber_NIST(train)));        
        
        temp11=polyfit(peak_positions_subpixel(i,test),wavenumber_NIST(test),6);
        temp12(i,:)=polyval(temp11,peak_positions_subpixel(i,train));
        error7(i,:) =mean(abs(temp12(i,:)-wavenumber_NIST(train)));        
        
        temp13=polyfit(peak_positions_subpixel(i,test),wavenumber_NIST(test),7);
        temp14(i,:)=polyval(temp13,peak_positions_subpixel(i,train));
        error8(i,:) =mean(abs(temp14(i,:)-wavenumber_NIST(train)));        
end

for i=1:100
std_error1(i)=std(wavenumber1(i,:)-wavenumber_NIST(train));      
std_orignial1(i)=std(temp2(i,:)-wavenumber_NIST(train));
std_orignial2(i)=std(temp4(i,:)-wavenumber_NIST(train));
std_orignial3(i)=std(temp6(i,:)-wavenumber_NIST(train));
std_orignial4(i)=std(temp8(i,:)-wavenumber_NIST(train));
std_orignial5(i)=std(temp10(i,:)-wavenumber_NIST(train));
std_orignial6(i)=std(temp12(i,:)-wavenumber_NIST(train));
std_orignial7(i)=std(temp14(i,:)-wavenumber_NIST(train));
end
STD_e1=mean(std_error1);STD_o1=mean(std_orignial1);STD_o2=mean(std_orignial2);STD_o3=mean(std_orignial3);
STD_o4=mean(std_orignial4);STD_o5=mean(std_orignial5);STD_o6=mean(std_orignial6);STD_o7=mean(std_orignial7);
for i=1:100
rmse_error1(i)=sqrt(mean((wavenumber1(i,:)-wavenumber_NIST(train)).^2));       
rmse_orignial1(i)=sqrt(mean((temp2(i,:)-wavenumber_NIST(train)).^2));
rmse_orignial2(i)=sqrt(mean((temp4(i,:)-wavenumber_NIST(train)).^2));
rmse_orignial3(i)=sqrt(mean((temp6(i,:)-wavenumber_NIST(train)).^2));
rmse_orignial4(i)=sqrt(mean((temp8(i,:)-wavenumber_NIST(train)).^2));
rmse_orignial5(i)=sqrt(mean((temp10(i,:)-wavenumber_NIST(train)).^2));
rmse_orignial6(i)=sqrt(mean((temp12(i,:)-wavenumber_NIST(train)).^2));
rmse_orignial7(i)=sqrt(mean((temp14(i,:)-wavenumber_NIST(train)).^2));
end
RMSE_e1=mean(rmse_error1);RMSE_o1=mean(rmse_orignial1);RMSE_o2=mean(rmse_orignial2);RMSE_o3=mean(rmse_orignial3);
RMSE_o4=mean(rmse_orignial4);RMSE_o5=mean(rmse_orignial5);RMSE_o6=mean(rmse_orignial6);RMSE_o7=mean(rmse_orignial7);
%% result=[Our method, 1st order, 2nd order,3rd order,4th order,5th order,6th order,7th order]
MAE_ace_g300_llh_peaks=[mean(error1_original1),mean(error2),mean(error3),mean(error4),mean(error5),mean(error6),mean(error7),mean(error8)];
RMSE_ace_g300_llh_peaks=[RMSE_e1,RMSE_o1,RMSE_o2,RMSE_o3,RMSE_o4,RMSE_o5,RMSE_o6,RMSE_o7];
STD_ace_g300_llh_peaks=[STD_e1,STD_o1,STD_o2,STD_o3,STD_o4,STD_o5,STD_o6,STD_o7];
