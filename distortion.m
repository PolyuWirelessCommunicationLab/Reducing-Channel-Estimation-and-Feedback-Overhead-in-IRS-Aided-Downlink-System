function Distortion = distortion(h,x)
% h is the vector(s) to be quantized
% x is the vector(s) in the codebook

% if size(h,2)<size(x,2)
%     Distortion=sum(abs(h).^2)-(abs(x'*h).^2).';
% elseif size(h,2)>size(x,2)
%     Distortion=sum(abs(h).^2,1)-abs(x'*h).^2;
% else
%     for i=1:size(h,2)
%         Distortion(i)=sum(abs(h(:,i)).^2)-abs(x(:,i)'*h(:,i)).^2;
%     end

Distortion=sum(abs(h-x).^2,1);

end