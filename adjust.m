%balance the weight between A and B
function [A2 B2]=adjust(A,B)
	wA=sum(sum(sum(A)))/numel(A);
	wB=sum(sum(sum(B)))/numel(B);
	A2=(A/wA)*sqrt(wA*wB);
	B2=(B/wB)*sqrt(wA*wB);
end