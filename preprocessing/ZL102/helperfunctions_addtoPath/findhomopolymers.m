function [number,h]=findhomopolymers(reads,runlength)
number=[];
Ind=[];
Ind5=[];

for i=1:size(reads,1)
Ind = [0, find(diff(reads(i,:))), length(reads(i,:))];
Len = diff(Ind);
Ind5 = Ind(Len >= runlength) + 1;
number(i)=size(Ind5,2);
end
h=hist(number,0:40);
end
