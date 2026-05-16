function [ Locat ] = FindLoca2( U1,Dataoedered )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
     LeftP=0;
     RirhtP=length(Dataoedered);
     while LeftP<RirhtP-1
           MidP=round((LeftP+RirhtP)/2);
           if U1<Dataoedered(MidP)
               RirhtP=MidP;
           else
               LeftP=MidP;
           end
           
     end
     Locat=LeftP;
     


end

