function [Tmin,Tsec] = sec2min(T)

Tmin = floor(T/60);
Tsec = round((T-Tmin*60)*10)/10;