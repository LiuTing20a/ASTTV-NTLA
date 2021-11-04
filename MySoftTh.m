%% This is soft thresholding operation
function X= MySoftTh(B,lambda)
X=sign(B).*max(0,abs(B)-(lambda));
end
