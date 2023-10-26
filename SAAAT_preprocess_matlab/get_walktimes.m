function[walktimes]=get_walktimes(velo,timewa)

 

SI=0.02; % sample interval; seconds
fc = 0.1;  %low pass filter cutoff frequency; Hz
fs=1/SI; %sampling frequency
start_cutoff=0.025; %velocity cutoff for walkbout start; m/s
end_cutoff=0.005; %velocity cutoff for walkbout end; m/s

 

walk_test=velo;%SessionWise(4).velocity;

 

%figure; plot(walk_test)

 

walk_rect=abs(walk_test);
%figure; plot(walk_rect)

 

[b,a] = butter(6,fc/(fs/2));

 

walk_smooth=filtfilt(b,a,walk_rect);
%
% figure; plot(walk_smooth)
% figure; histogram(walk_smooth)

 

%max(walk_smooth)

 

walk_starts=[];
for i=1:length(walk_smooth)-1
    if walk_smooth(i+1)>=start_cutoff && walk_smooth(i)<start_cutoff
        walk_starts=[walk_starts,i];
    end
end

 

walk_ends_all=[];
for i=1:length(walk_smooth)-1
    if walk_smooth(i+1)<=end_cutoff && walk_smooth(i)>end_cutoff
        walk_ends_all=[walk_ends_all,i];
    end
end
walk_ends_all=[walk_ends_all length(walk_smooth)+1];

 

walk_ends=nan(size(walk_starts));
for i=1:length(walk_starts)
    walk_ends_temp=walk_ends_all-walk_starts(i);
    walk_ends_temp(walk_ends_temp<=0)=[];
    walk_ends(i)=min(walk_ends_temp)+walk_starts(i);
end

 

remove_starts=[];
while length(remove_starts)>=1
    remove_starts=[];
    for i=2:length(walk_starts)
        if walk_starts(i)<=walk_ends(i-1)
            remove_starts=[remove_starts,i];
        end
    end
    walk_starts(remove_starts)=[];
    walk_ends(remove_starts)=[];
end
if walk_ends(length(walk_ends))>length(timewa)
    walk_ends(length(walk_ends))=length(timewa);
end

 

walktimes(:,1)=timewa(walk_starts);
walktimes(:,2)=timewa(walk_ends);

 

qq= unique(walktimes(:,2));

 

stored=[];
for kk=1:length(qq)
    if length(find(walktimes(:,2)==qq(kk)))>1
        a1=find(walktimes(:,2)==qq(kk));
        stored=[stored;a1(2:end)];
    end
end
walktimes(stored,:)=[];
end