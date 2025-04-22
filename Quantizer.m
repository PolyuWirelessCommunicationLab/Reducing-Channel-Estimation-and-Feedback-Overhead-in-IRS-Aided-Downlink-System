function vec_quan = Quantizer(vec_set, codebook)
vec_num=size(vec_set,2);
cb_size=size(codebook,2);
Distortion=zeros(cb_size,vec_num);
for i=1:cb_size
    Distortion(i,:)=sum(abs(vec_set-codebook(:,i)).^2,1);    
end
[~,quan_idx]=min(Distortion);
vec_quan=codebook(:,quan_idx);