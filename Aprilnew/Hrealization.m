L = 1000;
HZ = zeros(64,64,64,L);
AT = zeros(64,50,L);
AR = zeros(64,50,L);
for i = 1:L
    [HZ(:,:,:,i),AT(:,:,i),AR(:,:,i)]=OMPHWB(64,64);
end