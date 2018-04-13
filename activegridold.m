function [App, Dt]=activegridold(App, Dt, s, eps, lev)
%This function spits out the Dt, App, and y1 (finest resolution) structures
%for the active grid of a function with certain threshold on the wavelets.

%THIS IS THE OLD CODE FOR ACTIVE GRID BEFORE THE RESTRUCTURING ON 30/31
%MARCH


j=log2(length(s));
I2=find(abs(Dt')>eps); %finds points on the grid which are currently above the threshold. counts down column usually, but we want it to count down row.

len=length(App(1,:))*2;
%Now we need to find grid points which are adjacent in both scale and
%space. This is step 1 in Rong

%ADJACENT IN SPACE

%define which pts are the first pt on each scale
for i=1:lev
    firstpt(i)=(i-1)*2^(j-1)+1;
end

%define which pts are the last pt on each scale
for i=1:lev
    lastpt(i)=(i-1)*2^(j-1)+2^(j-i);
end

I3=I2;%store values in I3
for k=1:len*lev
     for i=1:lev
        if k>(len/2)*(i-1) && k<=(len/2)*(i) %determine which line the point is on.
            row=i;
        end
    end
    if ismember(k,I2)==1 && ismember(k,firstpt)==0 && ismember(k,lastpt)==0 %if k is significant but is not a first or last pt at each scale
        I3(length(I3)+1)=k-1; %add into the I2 array the wavelet prior in space
        I3(length(I3)+1)=k+1; %add into the I2 array the wavelet after in space
    end
    if ismember(k,I2)==1 && ismember(k,firstpt)==1 %if is the first point on some scale
        I3(length(I3)+1)=k+1; %add the second point on that scale
        I3(length(I3)+1)=k+2^(j-row)-1; %add the last point on that scale
    end
    if ismember(k,I2)==1 && ismember(k,lastpt)==1
        I3(length(I3)+1)=k-1;
        I3(length(I3)+1)=k-2^(j-row)+1;
    end
end


%ADJACENT IN SCALE
%We want to add indices to the I3 matrix which correlate to be adjacent in
%scale.
for k=1:len*lev
if ismember(k,I2)==1 %if we are in the I2 matrix
    for i=1:lev
        if k>(len/2)*(i-1) && k<=(len/2)*(i) %determine which line the point is on.
            row=i;
        end
    end
        m=k-(len/2)*(row-1); %gives position on line
        if mod(m,2)==0 && row~=1 && row~=lev %if an even point on the line, and the line isn't the first or last
            I3(length(I3)+1)=k+(len/2-m) + m/2; %goes down one line
            I3(length(I3)+1)=k-m-(len/2)+2*m; % goes up one line
        end
        if mod(m,2)~=0 && row~=1 && row~=lev || row==lev %if odd on a middle line, or if on the last line
            I3(length(I3)+1)=k-m-(len/2)+2*m; % goes up one line
        end
        if row==1 && mod(m,2)==0
             I3(length(I3)+1)=k+(len/2-m) + m/2; %goes down one line
        end            
    end
end
%Remember, we have to transpose it because we counted the transposed way.    
B=Dt';
I4=0;
for v=1:lev*(len/2)
    if ismember(v,I3)==0
        I4(length(I4)+1)=v; %store all wavelet values not in safety zone
    end
end
I4=I4(2:end)';
B(I4)=zeros(size(I4));
Dt=B';


%MODIFY THE APP MATRIX
%find wavelets at coarsest level
for i=1:(len/(2^lev))
I5(i)=(lev-1)*(len/2)+i;
end

I6=0; %initialize the list that we will store the values at the finest resolution level that we are going to keep

%find those in perfect reconstruction zone above each active wavelet in I2
for k=1:len*lev
if ismember(k,I2)==1 %if we are in the I2 matrix
    for i=1:lev
        if k>(len/2)*(i-1) && k<=(len/2)*(i) %determine which line the point is on.
            row=i;
        end
    end
        m=k-(len/2)*(row-1); %gives position on line
        if row~=1 && m~=len/2^row %if not the last wavelet on a line and not on the first line
            I5(length(I5)+1)=k-(len/2) +m; %gets the scaling function to the left of active wavelet
            I5(length(I5)+1)=k-(len/2)+m+1;
        end
        if row~=1 && m==len/2^row %if not on the first line, but are the last wavelet on a line
            I5(length(I5)+1)=k-(len/2)+m;
            I5(length(I5)+1)=k-(len/2)-m+1;
        end  
%         if row==1 && m==len/2 %if the last wavelet at the first level of decomposition, we take the first and last pt at highest level of decomp
%             I6(length(I6)+1)=len;
%             I6(length(I6)+1)=1; WE DON'T NEED TO DO THIS STEP. SEE THE
%             STAR NOTE ON THE BACK OF THE 25/26 MARCH PAGE
%         end
%         if row==1 && m~=len/2 %if on the first level of decomposition and not a centre wavelet
%             I6(length(I6)+1)=2*m-1;
%             I6(length(I6)+1)=2*m+1;
%         end
    end
end
%Remember, we have to transpose it because we counted the transposed way. 
%Make everything in the App matrix 0 unless it is in the perfect
%reconstruction zone or the coarsest level
C=App';
I7=0;
for v=1:lev*(len/2)
    if ismember(v,I5)==0
        I7(length(I7)+1)=v; %store all scaling function values not in perfect reconstruction zone or in coarsest level
    end
end
I7=I7(2:end)';
C(I7)=zeros(size(I7));
App=C';

%make everything in the finest level 0 unless it is in the perfect
%reconstruction zone
% I8=0; %SEE THE STAR ON THE BACK OF THE 25/26 MARCH PAPER
% for v=1:len
%     if ismember(v,I6)==0
%         I8(length(I8)+1)=v;
%     end
% end
% I8=I8(2:end);
% y1=s; %for now, we aren't going to throw away the whole solution
% y1(I8)=zeros(size(I8));
end
