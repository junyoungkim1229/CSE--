clc;
clear;

A=randi([1,100],1,100)
%randi:지정한 범위 내의 임의의 난수 생성
%randi([imin,imax],m(행),n(열));

B=sort(A)
%sort:오름차순 정렬,
%내림차순은 B=sort(A,'descend')


