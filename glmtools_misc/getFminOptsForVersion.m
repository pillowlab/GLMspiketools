function opts = getFminOptsForVersion(v)
% opts = getFminOptsForVersion(v)
%
% Gets the options for fminunc needed for user's version of MATLAB
%
% Note of complaint: it's stupid that I had to write this function, but
% matlab has incompatible syntax for fminunc options for different versions
% (earlier versions of matlab will crash with options required by later
% versions).  The cutoff point for my two versions was 9.1 vs. 9.3, but if
% you're getting warnings or crashes, please report to pillow@princeton.edu
% or make an issue report on github

if str2num(v(1:3))>9.2 % check if version is before V 9.2
    opts = {'algorithm','trust-region','Gradobj','on','Hessian','on'};
else
    opts = {'Gradobj','on','Hessian','on'};
end
