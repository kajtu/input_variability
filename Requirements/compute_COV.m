function [ COV ] = compute_COV( x,y )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%COV: Within-subject coefficient of variation (%)
%x: Dataset 1 (data from analysis/observer1)
%y: Dataset 2 (data from analysis/observer2)

s2=(x-y).^2/2;
m=(x+y)./2;
s2m2=s2./m.^2;
means2m2=mean(s2m2);
COV=sqrt(means2m2)*100;

end

