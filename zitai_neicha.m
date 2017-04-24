zitai=dlmread('./time_rotate.txt');
load('NFB_line.mat');
f_end=find(B_Line(:,1)==F_Line(end,1));
time=[F_Line(:,8);B_Line(f_end+1:end,8)];

for i=1:length(time)
[array,index]=sort(abs(zitai(:,1)-time(i)));
fai(i)=lagrange(zitai(index(1:2),1),zitai(index(1:2),2),time(i));
omg(i)=lagrange(zitai(index(1:2),1),zitai(index(1:2),3),time(i));
kappa(i)=lagrange(zitai(index(1:2),1),zitai(index(1:2),4),time(i));
end

time_rotate_chazhi(:,1)=time-time(1);
time_rotate_chazhi(:,2)=fai;
time_rotate_chazhi(:,3)=omg;
time_rotate_chazhi(:,4)=kappa;
fileID=fopen('time_rotate_chazhi.txt','w');
for i=1:length(time)
    fprintf(fileID,'%.8f %.20f %.20f %.20f\n',time_rotate_chazhi(i,1),...
        time_rotate_chazhi(i,2),time_rotate_chazhi(i,3),time_rotate_chazhi(i,4));
end
fclose(fileID);
