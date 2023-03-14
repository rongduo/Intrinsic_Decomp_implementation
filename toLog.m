%convert the A to B=log(A+epsilon) robustly
function B=toLog(A)
	epsilon=max(max(max(A)))/1000;
	B=log(A+epsilon);
end