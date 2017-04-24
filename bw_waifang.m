function [y,W,Xs,Ys,Zs]=bw_waifang(u,v)
load('NFB_line.mat');
NC2ECF=dlmread('time_rotate_chazhi.txt');
OXYZ=dlmread('TH01-01_P201305251109298_1A_SXZ.JDG');
first_line_num=B_Line(1,1);
line_start_time=154851190.4578571;
y0=6002.34;
pixel_size=6.5*1e-6;
y=(y0-v)*pixel_size;
hang_num=floor(u);
u_hang=[hang_num,hang_num+1]';
hang_index=hang_num-first_line_num+1;
hang_time=[B_Line(hang_index,8),B_Line(hang_index+1,8)];
line_time=interp1(u_hang,hang_time,u);
line=round(((line_time-line_start_time)/(NC2ECF(end,1)/(length(NC2ECF)-1))))+1;
relevent_time=[NC2ECF(line-1,1),NC2ECF(line,1),NC2ECF(line+1,1),NC2ECF(line+2,1)]';
fi0=[NC2ECF(line-1,2),NC2ECF(line,2),NC2ECF(line+1,2),NC2ECF(line+2,2)]';
si0=[NC2ECF(line-1,3),NC2ECF(line,3),NC2ECF(line+1,3),NC2ECF(line+2,3)]';
ka0=[NC2ECF(line-1,4),NC2ECF(line,4),NC2ECF(line+1,4),NC2ECF(line+2,4)]';
Xs_0=[OXYZ(line-1,2),OXYZ(line,2),OXYZ(line+1,2),OXYZ(line+2,2)]';
Ys_0=[OXYZ(line-1,3),OXYZ(line,3),OXYZ(line+1,3),OXYZ(line+2,3)]';
Zs_0=[OXYZ(line-1,4),OXYZ(line,4),OXYZ(line+1,4),OXYZ(line+2,4)]';
fi=lagrange(relevent_time,fi0,line_time-line_start_time);
si=lagrange(relevent_time,si0,line_time-line_start_time);
ka=lagrange(relevent_time,ka0,line_time-line_start_time);
W1=Ro1(fi,si,ka);
Xs=lagrange(relevent_time,Xs_0,line_time-line_start_time);
Ys=lagrange(relevent_time,Ys_0,line_time-line_start_time);
Zs=lagrange(relevent_time,Zs_0,line_time-line_start_time);
R3=[0.9061989451  -0.0000139277  0.4228515957
    0.0000250280  0.9999999995  -0.0000206472 
    -0.4228515957  0.0000293041  0.9061989451];
%W=R3'*W1;
W=W1*R3';
end
