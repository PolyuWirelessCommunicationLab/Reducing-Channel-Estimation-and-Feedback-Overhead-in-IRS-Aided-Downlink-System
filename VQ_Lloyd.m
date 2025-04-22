function [codebook, Distortion_all] = VQ_Lloyd(training_set, bit)
% Lloyd algorithm to generate the quantization codebook
vec_len=size(training_set,1);
vec_num=size(training_set,2);
size_cb=2^bit;
%initialize
ini_idx=randperm(vec_num,size_cb);
ini_cb=training_set(:,ini_idx);
codebook=ini_cb;
Distortion=0;
Distortion_all=[];
max_iter=5e2;
epsilon=1e-3;
iter=0;
while (iter<max_iter)
    Reg_div=cell(1,size_cb);
    for vec_idx=1:vec_num
        dist=distortion(training_set(:,vec_idx),codebook);
        [~,reg_idx]=min(dist);
        Reg_div{reg_idx}=[Reg_div{reg_idx}, vec_idx];            
    end
    Dist=zeros(1,size_cb);
    for cb_idx=1:size_cb    
        %For MSE distortion
        codebook(:,cb_idx)=mean(training_set(:,Reg_div{cb_idx}),2);
        Dist(cb_idx)=sum(distortion(training_set(:,Reg_div{cb_idx}),codebook(:,cb_idx)));%/size(Reg_div{cb_idx},2);
    end
    Distortion_new=sum(Dist)/vec_num;
    if abs((Distortion_new-Distortion)/Distortion)<epsilon
        Distortion_all=[Distortion_all,Distortion_new];
        break;
    end
    Distortion=Distortion_new;
    iter=iter+1;
    Distortion_all=[Distortion_all,Distortion];
end

end
