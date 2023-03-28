function chosen_k=k_chosen(peak_positions_pixel,peak_numbers,rawdata)
wind1=[3,3,3,3,3,3,3,2,2,3,3,3,3,3,3,3,3,3,2,3]; % The left part of window
wind2=[3,3,3,3,3,3,3,2,2,3,3,3,3,3,3,3,3,3,3,3]; % The right part of window
p=peak_positions_pixel-wind1;
q=peak_positions_pixel+wind2; 
Min1=zeros(100,peak_numbers);
chosen_k=zeros(100,peak_numbers);
chosen_P=zeros(100,peak_numbers);
peak_positions_subpixel=zeros(100,peak_numbers);
parfor i=1:100
    for j=1:peak_numbers
        Min1(i,j)=1e5; %check
        for k=0.1:0.1:100 
        y=rawdata(i,p(i,j):q(i,j));   
        x=p(i,j):1:q(i,j);
        P2=peak_positions_pixel(i,j);P3=k;C=0;
        [height,pos]= max(y);
        P1=P3*height;
        P0=[P1,P2,P3,C];
        [yprime2,P,resnorm2,residual2] = lorentzfit(x,y,P0);
        [peak_positions_subpixel(i,j)]=P(2);
        xx=p(i,j):0.01:q(i,j);
        YPRIME = P(1)./((xx - P(2)).^2 + P(3)) + P(4);
        line=[xx',YPRIME'];
        [index,temp1]=find(x==line(:,1));
        test=sum(abs(line(index,2)-y'));
            if test<=Min1(i,j)
            Min1(i,j)=test;
            chosen_k(i,j)=k;
            chosen_P(i,j)=mean(abs(P));
            end
        end
    end
end

end