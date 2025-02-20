%Fan-Beam ART Siddon V7
%----------------------------------Description---------------------------------%
%系统矩阵:Siddon(预先计算系统矩阵,占用内存大,速度快)
%重建算法:ART
%基础程序:FBAS_V6
%------------------------------------------------------------------------------%
tic  %开始计时
%----------------------------------Parameters----------------------------------%
det   = 'A';   %探测器排列方式(A=等角;L=等距)
pix   = 512;   %图像尺寸(单边像素数量)
sod   = 726;   %射线源-旋转中心距离(sod>0.5*pix/(sin(fan/2)*sin(45°)))
fan   = 60;    %扇形束张角(角度制)
anum  = 1200;  %旋转角度数量(射线源沿逆时针方向旋转)
dnum  = 801;   %探测器数量(探测器沿逆时针方向编号)
fidx  = 1;     %滤波器编号(1=Ramp;2=Hanning;3=Hamming;4=Cosine;5=Triang;6=BartHann)
inum  = 6;     %迭代次数
Relax = 0.1;   %松弛因子(控制迭代的步长)
TV    = 0;     %TV正则(1=使用;0=不使用)
Mode  = 1;     %迭代重建模式(0=基于空图;1=基于FBP)
Range = 1;     %0=每轮约束<0的值;1=每轮约束<0和>1的值;其他=无操作
%----------------------------------File Name-----------------------------------%
FN = strcat('FBAS_V7',det,'p',string(pix),'s',string(sod),'f',string(fan*10),'a',string(anum),...
            'd',string(dnum),'f',string(fidx),'i',string(inum),'r',string(Relax));
%------------------------------------------------------------------------------%
img   = phantom(pix);  %重建对象
iseq  = img(:);  %重建对象(列向量)
hpix  = 0.5*pix;
FBP   = zeros(pix*pix,1);  %FBP数据
mask  = zeros(pix*pix,1);  %背景区域(背景区域>0,目标区域=0)
if strcmp(det,'A')  %射线的初始角度(等角排布)
    iba = transpose(linspace(-fan/2,fan/2,dnum));  
elseif strcmp(det,'L')  %射线的初始角度(等距排布)
    iba = transpose(atand(linspace(tand(-fan/2),tand(fan/2),dnum)));
end
SIDi  = zeros(pix*2-1,dnum*anum);  %Siddon index(射线穿过像素的编号)
SIDw  = zeros(pix*2-1,dnum*anum);  %Siddon weight(射线穿过像素的权重)
SIDn  = zeros(dnum,anum);  %Siddon number(射线穿过像素的数量)
sino  = zeros(dnum,anum);  %Sinogram
%-------------------------------Weight and Filter------------------------------%
wseq  = 1./cosd(iba);  %加权序列
flen  = 16384;  %滤波器长度
fseq  = linspace(-1,1,flen);
win   = ones(10,flen);  %窗口函数
win(2,:) = hann(flen);
win(3,:) = hamming(flen);
win(4,:) = cos(pi*fseq/2);
win(5,:) = triang(flen);
win(6,:) = barthannwin(flen);
filt  = transpose(abs(fseq).*win(fidx,:));
%----------------------------Calculate System Matrix---------------------------%
for a = 1:anum
    sa = (a-1)*360/anum;  %source angle(角度制)
    ba = iba+sa;  %beam angle(角度制)
    x0 = sod*cosd(sa);  %射线源的X轴坐标
    y0 = sod*sind(sa);  %射线源的Y轴坐标
    for d = 1:dnum
        k = tand(ba(d));  %射线的斜率
        %--------射线垂直于X轴--------%
        if abs(k) == inf
            if abs(x0) <= hpix  %判断射线是否穿过图像
                col = ceil(x0)+hpix;  %确定射线对应的列
                SIDi(1:pix,(a-1)*dnum+d) = (col-1)*pix+1:1:col*pix;
                SIDw(1:pix,(a-1)*dnum+d) = ones(pix,1);
                SIDn(d,a) = pix;
            end
        %--------射线垂直于Y轴--------%
        elseif k == 0
            if abs(y0) <= hpix  %判断射线是否穿过图像
                row = -floor(y0)+hpix;  %确定射线对应的行
                SIDi(1:pix,(a-1)*dnum+d) = row:pix:(pix-1)*pix+row;
                SIDw(1:pix,(a-1)*dnum+d) = ones(pix,1);
                SIDn(d,a) = pix;
            end
        %--------射线处于其他角度--------%
        else
            xpart = x0-y0/k;  %x=y/k+xpart(已知y,求x)
            ypart = y0-x0*k;  %y=x*k+ypart(已知x,求y)
            xT = hpix/k+xpart;  %x top
            xB = -hpix/k+xpart;  %x bottom
            passR = xT>hpix && xB>hpix;  %判断射线是否从图像右侧经过
            passL = xT<-hpix && xB<-hpix;  %判断射线是否从图像左侧经过
            if passR+passL ~= 0  %射线与图像无交点时跳过本次循环
                continue;
            end
            COS = abs(cosd(ba(d)));  %射线角度余弦值
            %--------k<0,射线方向右下/左上--------%
            if k < 0
                xL = max(-hpix,xT);  %射线与图像左交点X坐标
                xR = min(hpix,xB);  %射线与图像右交点X坐标
            %--------k>0,射线方向右上/左下--------%
            else
                xL = max(-hpix,xB);  %射线与图像左交点X坐标
                xR = min(hpix,xT);  %射线与图像右交点X坐标
            end
            yLR = [xL*k+ypart,xR*k+ypart];  %射线与图像左右交点Y坐标
            ymid = floor(min(yLR))+1:1:ceil(max(yLR))-1;  %左右交点之间的整数Y坐标
            xcord = unique([xL,floor(xL)+1:1:ceil(xR)-1,xR,ymid/k+xpart]);  %所有交点的X坐标
            ycord = xcord*k+ypart;  %所有交点的Y坐标
            xnum = length(xcord);
            rnow = hpix-floor((ycord(1:xnum-1)+ycord(2:xnum))*0.5);  %射线穿过图像的行
            cnow = hpix+1+floor(xcord(1:xnum-1));  %射线穿过图像的列(xcord从小到大,简化计算)
            SIDi(1:xnum-1,(a-1)*dnum+d) = (cnow-1)*pix+rnow;
            SIDw(1:xnum-1,(a-1)*dnum+d) = transpose((xcord(2:xnum)-xcord(1:xnum-1))/COS);
            SIDn(d,a) = xnum-1;
        end
    end
end
%-------------------------------FBP at each angle------------------------------%
for a = 1:anum
    %--------Projection--------%
    PD = zeros(dnum,1);  %投影数据
    for d = 1:dnum
        index = SIDi(1:SIDn(d,a),(a-1)*dnum+d);  %射线穿过像素的编号
        PD(d) = sum(iseq(index).*SIDw(1:SIDn(d,a),(a-1)*dnum+d));  %记录射线的投影值
        if PD(d) == 0  %投影值=0,记录射线穿过像素的编号
            mask(index) = 1;
        end
    end
    %--------Weighting/Filtering/Backprojection--------%
    wPD   = PD;  %Siddon算法不加权的重建效果更好
    fftp  = fftshift(fft(wPD,flen));  %fftshift将fft变为负高频-低频-高频形式
    ifftp = real(ifft(ifftshift(fftp.*filt)));  %滤波+傅里叶反变换
    fwPD  = ifftp(1:dnum);
    for d = 1:dnum  %反投影
        index = SIDi(1:SIDn(d,a),(a-1)*dnum+d);  %射线穿过像素的编号
        FBP(index) = FBP(index)+fwPD(d)*SIDw(1:SIDn(d,a),(a-1)*dnum+d);
    end
    sino(:,a) = PD;  %将当前角度的投影数据记录到Sinogram中
end
%---------------------------Iterative Reconstruction---------------------------%
FBPf = min(max(FBP*pi/(2*anum),0),1);  %将FBP重建数据规格化并约束在[0,1]的范围
FBPo = FBPf.*(1-logical(mask));  %FBP重建数据(去除背景信息)
if Mode == 0
    ARTf = zeros(pix*pix,1);  %基于空图重建
elseif Mode == 1
    ARTf = FBPo;  %基于FBP结果重建
end
for i = 1:inum
    for a = 1:anum
        for d = 1:dnum
            index = SIDi(1:SIDn(d,a),(a-1)*dnum+d);  %射线穿过像素的编号
            if sino(d,a) == 0  %原始投影值=0,射线穿过像素赋值0
                ARTf(index) = 0;
            else
                weight = SIDw(1:SIDn(d,a),(a-1)*dnum+d);  %当前射线穿过像素的权重
                PD = sum(ARTf(index).*weight);  %当前射线的投影值
                ARTf(index) = ARTf(index)+Relax*weight*(sino(d,a)-PD)/sum(weight.^2);  %更新图像
            end
        end
    end
    if Range == 0
        ARTf = max(ARTf,0);  %约束<0的值
    elseif Range == 1
        ARTf = min(max(ARTf,0),1);  %约束<0和>1的值
    end
    %--------TV Regularization--------%
    if TV == 1
        imgTV = reshape(ARTf,pix,pix);
        [gradx,grady] = gradient(imgTV);
        gmag = sqrt(gradx.^2+grady.^2+1e-8);  %gradient magnitude
        div = divergence(gradx./gmag,grady./gmag);
        imgTV = imgTV+0.0025*div;
        ARTf = imgTV(:);
    end
end
%--------------------------------Reconstruction--------------------------------%
ARTf = min(max(ARTf,0),1);  %将ART重建数据约束在[0,1]的范围
ARTo = ARTf.*(1-logical(mask));  %ART重建数据(去除背景信息)
IMGf = reshape(ARTf,pix,pix);  %ART重建图像(完整)
IMGo = reshape(ARTo,pix,pix);  %ART重建图像(去除背景信息)
%------------------------Evaluate Reconstruction Quality-----------------------%
S_MSE_f  = mean((ARTf-img(:)).^2);
S_MSE_o  = mean((ARTo-img(:)).^2);
S_PSNR_f = psnr(IMGf,img);
S_PSNR_o = psnr(IMGo,img);
S_SSIM_f = ssim(IMGf,img);
S_SSIM_o = ssim(IMGo,img);
Sum = [S_MSE_f,S_PSNR_f,S_SSIM_f,S_MSE_o,S_PSNR_o,S_SSIM_o];
%------------------------------------------------------------------------------%
%subplot(1,2,1),imshow(img),title('Original')
%subplot(1,2,2),imshow(IMGf),title('ART')
pseq = linspace(1,pix,pix);
hsec = pix/2;  %horizontal section
vsec = pix/2;  %vertical section
subplot(2,1,1),plot(pseq,img(hsec,:),pseq,IMGf(hsec,:),'r');
subplot(2,1,2),plot(pseq,img(:,vsec),pseq,IMGf(:,vsec),'r');
%imwrite(IMGf,strcat(FN,'.bmp'));
%imwrite(IMGo,strcat(FN,'NBG.bmp'));
%------------------------------------------------------------------------------%
toc  %结束计时
%------------------------------------------------------------------------------%