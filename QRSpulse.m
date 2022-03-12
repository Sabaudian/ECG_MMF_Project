function [QRS_wave]=QRSpulse(qrswidth,bpm,fs,amp)
%QRSpulse(qrswidth,bpm,fs,amp) artifical QRS pulse generator
%   [QRS_wave]=QRSPULSE(fs,bpm,amp,qrswidth) returns a waveform vector.
%   fs is the sample frequency, 100 to 500 Hz typical.
%   bpm is the number of beats per minute.
%   amp is the amplitude of the QRS in micro volts, 500 to 5000 uV.
%   qrswidth in seconds, the time period between QRS.

%creates an artifical ECG signal



if (exist('fs') ~= 1)  fs=  200;   end  %default value, Hz
if (exist('bpm') ~= 1)  bpm =  72;   end %default value, beats per minute
if (exist('amp') ~= 1)  amp = 1000;   end %default value, micro volts
if (exist('qrswidth') ~= 1)  
    d = 0.135;   %.07 to .120 seconds, QRS width
else
    d=qrswidth; %.07 to .120 seconds, QRS width
end %default value, micro volts


%things that can change
%d=0.070; %.07 to .120 seconds, QRS width
at=500; %amplitude of t-wave, 400 to 1200 uv

RR=(60/bpm); %RR interval
d1=0.4375*d;
d2=0.5*d;
d3=d-(d1+d2);
dt=0.180; %seconds
qt=0.35;
t_line1=0:1/fs:RR+1; QRS_wave=zeros( size(t_line1) ); 
deadspace=RR-qt; 
if deadspace < 0 
    err_msg=['Bpm must be less than ' int2str(60/qt) '. '];
    error(err_msg); 
end


%Segment 1 bottom-top (Q-R)
t1=0; t2=t1+d1;
i_t1=time2index(t1,t_line1); i_t2=time2index(t2,t_line1);
bottom=0; top=0.875*amp;
QRS1=bottom:(top-bottom)/(i_t2-i_t1):top;
QRS_wave(i_t1:i_t2)=QRS1;

%Segment 2 top-bottom (R-?)
t1=t2; t2=t1+d2;
i_t1=time2index(t1,t_line1); i_t2=time2index(t2,t_line1);
bottom=-.125*amp;
QRS1=top:(bottom-top)/(i_t2-i_t1):bottom;
QRS_wave(i_t1:i_t2)=QRS1;


%Segment 3 bottom_top (?-S) 
t1=t2; t2=t1+d3; 
i_t1=time2index(t1,t_line1); i_t2=time2index(t2,t_line1);
top=0;
if (i_t2-i_t1 >0)
    QRS1=bottom:(top-bottom)/(i_t2-i_t1):top;
    QRS_wave(i_t1:i_t2)=QRS1;
end

%Segment 4 horizontal line
t1=t2; t2=t1+qt-(dt+t2);
i_t1=time2index(t1,t_line1); i_t2=time2index(t2,t_line1);
bottom=0;

%Segment 5, half circle
t1=t2; t2=t1+dt;
i_t1=time2index(t1,t_line1); i_t2=time2index(t2,t_line1);
t=-1:2/(i_t2-i_t1):1;
QRS1=at*sqrt(1-(t).^2);
QRS_wave(i_t1:i_t2)=QRS1;

%Segment 6, rest of deadspace
t1=t2; t2=t1+deadspace;
i_t1=time2index(t1,t_line1); i_t2=time2index(t2,t_line1);
QRS_wave=QRS_wave(1:i_t2);


%stem(QRS_wave); %one cycle


function index=time2index(t, t_line1)
%TIME2INDEX converts time (s) to an index value

indexArray=find(t_line1>=t);
index=indexArray(1); 

