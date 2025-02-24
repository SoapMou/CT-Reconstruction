%Fan-Beam ART Siddon V6
%----------------------------------Description---------------------------------%
%系统矩阵:Siddon(即时计算系统矩阵,占用内存小,速度慢)
%重建算法:ART
%基础程序:FBARTV5
%------------------------------------------------------------------------------%
tic  %开始计时
%----------------------------------Parameters----------------------------------%
Gdet  = 'A';   %探测器排列方式(A=等角;L=等距)
Gpix  = 512;   %图像尺寸(单边像素数量)
Gsod  = 726;   %射线源-旋转中心距离(sod>0.5*pix/(sin(fan/2)*sin(45°)))
Gfan  = 60;    %扇形束张角(角度制)
Ganum = 600;  %旋转角度数量(射线源沿逆时针方向旋转)
Gdnum = 801;   %探测器数量(探测器沿逆时针方向编号)
Gfidx = 1;     %滤波器编号(1=Ramp;2=Hanning;3=Hamming;4=Cosine)
Ginum = 1;     %迭代次数
Relax = 0.1;   %松弛因子(控制迭代的步长)
TV    = 0;     %TV正则(1=使用;0=不使用)
Mode  = 1;     %迭代重建模式(0=基于空图;1=基于FBP)
Range = 1;     %0=每轮约束<0的值;1=每轮约束<0和>1的值;其他=无操作
Gimg  = phantom(Gpix);  %重建对象
%-----------------------------------File Name----------------------------------%
FN = strcat('FBAS_V6',Gdet,'p',string(Gpix),'s',string(Gsod),'f',string(Gfan*10),'a',string(Ganum),...
            'd',string(Gdnum),'f',string(Gfidx),'i',string(Ginum),'r',string(Relax));
%------------------------------------------------------------------------------%
GFBP  = zeros(Gpix*Gpix,1);  %FBP数据
Gmask = zeros(Gpix*Gpix,1);  %背景区域(背景区域>0,目标区域=0)
if strcmp(Gdet,'A')  %射线的初始角度(等角排布)
    Giba = transpose(linspace(-Gfan/2,Gfan/2,Gdnum));  
elseif strcmp(Gdet,'L')  %射线的初始角度(等距排布)
    Giba = transpose(atand(linspace(tand(-Gfan/2),tand(Gfan/2),Gdnum)));
end
sino  = zeros(Gdnum,Ganum);  %Sinogram
%------------------------------------Filter------------------------------------%
Gflen = 4096;  %滤波器长度
win   = ones(10,Gflen);  %窗口函数
win(2,:) = hann(Gflen);
win(3,:) = hamming(Gflen);
win(4,:) = cos(pi*linspace(-1,1,Gflen)/2);
Gfilt = transpose(abs(linspace(-1,1,Gflen)).*win(Gfidx,:));
%-------------------------------FBP at each angle------------------------------%
parfor a = 1:Ganum
    %--------全局变量->局部变量--------%
    pix  = Gpix;
    dnum = Gdnum;
    img  = Gimg;
    iseq = img(:);
    sa   = (a-1)*360/Ganum;  %source angle(角度制)
    ba   = Giba+sa;  %beam angle(角度制)
    x0   = Gsod*cosd(sa);  %射线源的X轴坐标
    y0   = Gsod*sind(sa);  %射线源的Y轴坐标
    %--------其他局部变量--------%
    hpix = pix/2; 
    FBP  = zeros(pix*pix,1);  %FBP数据
    PD   = zeros(dnum,1);  %投影数据
    mask = zeros(pix*pix,1);  %背景区域(背景区域=1,目标区域=0)
    SIDi = zeros(pix*2-1,dnum);  %Siddon index(射线穿过像素的编号)
    SIDw = zeros(pix*2-1,dnum);  %Siddon weight(射线穿过像素的权重)
    SIDn = zeros(dnum,1);  %Siddon number(射线穿过像素的数量)
    %--------Projection of each beam--------%
    for d = 1:dnum
        k = tand(ba(d));  %射线的斜率
        %--------射线垂直于X轴--------%
        if abs(k) == inf
            if abs(x0) <= hpix  %判断射线是否穿过图像
                col = ceil(x0)+hpix;  %确定射线对应的列
                PD(d) = sum(img(:,col));  %记录射线的投影值
                SIDi(1:pix,d) = (col-1)*pix+1:1:col*pix;  %记录射线穿过像素的编号
                SIDw(1:pix,d) = ones(pix,1);  %记录射线穿过像素的权重
                SIDn(d) = pix;  %射线穿过了pix个像素
                if PD(d) == 0  %投影值=0,记录射线穿过像素的编号
                    mask(SIDi(1:pix,d)) = 1;
                end
            end
        %--------射线垂直于Y轴--------%
        elseif k == 0
            if abs(y0) <= hpix  %判断射线是否穿过图像
                row = -floor(y0)+hpix;  %确定射线对应的行
                PD(d) = sum(img(row,:));  %记录射线的投影值
                SIDi(1:pix,d) = row:pix:(pix-1)*pix+row;  %记录射线穿过像素的编号
                SIDw(1:pix,d) = ones(pix,1);  %记录射线穿过像素的权重
                SIDn(d) = pix;  %射线穿过了pix个像素
                if PD(d) == 0  %投影值=0,记录射线穿过像素的编号
                    mask(SIDi(1:pix,d)) = 1;
                end
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
            index = (cnow-1)*pix+rnow;  %射线穿过像素的编号
            weight = transpose((xcord(2:xnum)-xcord(1:xnum-1))/COS);  %射线穿过像素的权重
            %--------计算投影值--------%
            PD(d) = sum(iseq(index).*weight);  %记录射线的投影值
            SIDi(1:xnum-1,d) = index;
            SIDw(1:xnum-1,d) = weight;
            SIDn(d) = xnum-1;  %射线穿过了xnum-1个像素
            if PD(d) == 0  %投影值=0,记录射线穿过像素的编号
                mask(index) = 1;
            end
        end
    end
    %--------Filtering/Backprojection--------%
    fftp  = fftshift(fft(PD,Gflen));  %fftshift将fft变为负高频-低频-高频形式
    ifftp = real(ifft(ifftshift(fftp.*Gfilt)));  %滤波+傅里叶反变换
    fPD   = ifftp(1:dnum);
    for d = 1:dnum  %反投影
        FBP(SIDi(1:SIDn(d),d)) = FBP(SIDi(1:SIDn(d),d))+fPD(d)*SIDw(1:SIDn(d),d);
    end
    GFBP = GFBP+FBP;  %将当前角度的FBP数据累加到GFBP中
    Gmask = Gmask+mask;  %将当前角度的mask累加到Gmask中
    sino(:,a) = PD;  %将当前角度的投影数据记录到Sinogram中
end
%---------------------------Iterative Reconstruction---------------------------%
FBPf = min(max(GFBP*pi/(2*Ganum),0),1);  %将FBP重建数据规格化并约束在[0,1]的范围
FBPo = FBPf.*(1-logical(Gmask));  %FBP重建数据(去除背景信息)
pix  = Gpix;
hpix = 0.5*pix;
if Mode == 0
    ARTf = zeros(pix*pix,1);  %基于空图重建
elseif Mode == 1
    ARTf = FBPo;  %基于FBP结果重建
end
for i = 1:Ginum
    for a = 1:Ganum
        sa = (a-1)*360/Ganum;  %source angle(角度制)
        ba = Giba+sa;  %beam angle(角度制)
        x0 = Gsod*cosd(sa);  %射线源的X轴坐标
        y0 = Gsod*sind(sa);  %射线源的Y轴坐标
        %--------Update image at each beam--------%
        for d = 1:Gdnum  %每分析一条射线更新一次图像
            k = tand(ba(d));  %射线的斜率
            %--------射线垂直于x轴--------%
            if abs(k) == inf
                if abs(x0) <= hpix  %判断射线是否穿过图像
                    col = ceil(x0)+hpix;  %确定射线对应的列
                    index = (col-1)*pix+1:1:col*pix;  %射线穿过像素的编号
                    PD = sum(ARTf(index));  %射线的投影值
                    ARTf(index) = ARTf(index)+Relax*(sino(d,a)-PD)/pix;  %更新图像
                end
            %--------射线垂直于y轴--------%
            elseif k == 0
                if abs(y0) <= hpix  %判断射线是否穿过图像
                    row = -floor(y0)+hpix;  %确定射线对应的行
                    index = row:pix:(pix-1)*pix+row;  %射线穿过像素的编号
                    PD = sum(ARTf(index));  %射线的投影值
                    ARTf(index) = ARTf(index)+Relax*(sino(d,a)-PD)/pix;  %更新图像
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
                index = (cnow-1)*pix+rnow;  %射线穿过像素的编号
                weight = transpose((xcord(2:xnum)-xcord(1:xnum-1))/COS);  %射线穿过像素的权重
                PD = sum(ARTf(index).*weight);  %射线的投影值
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
        imgTV = reshape(ARTf,Gpix,Gpix);
        [gradx,grady] = gradient(imgTV);
        gmag = sqrt(gradx.^2+grady.^2+1e-8);  %gradient magnitude
        div = divergence(gradx./gmag,grady./gmag);
        imgTV = imgTV+0.0025*div;
        ARTf = imgTV(:);
    end
end
%--------------------------------Reconstruction--------------------------------%
ARTf = min(max(ARTf,0),1);  %将ART重建数据约束在[0,1]的范围
ARTo = ARTf.*(1-logical(Gmask));  %ART重建数据(去除背景信息)
IMGf = reshape(ARTf,Gpix,Gpix);  %ART重建图像(完整)
IMGo = reshape(ARTo,Gpix,Gpix);  %ART重建图像(去除背景信息)
%------------------------Evaluate Reconstruction Quality-----------------------%
S_MSE_f  = mean((ARTf-Gimg(:)).^2);
S_MSE_o  = mean((ARTo-Gimg(:)).^2);
S_PSNR_f = psnr(IMGf,Gimg);
S_PSNR_o = psnr(IMGo,Gimg);
S_SSIM_f = ssim(IMGf,Gimg);
S_SSIM_o = ssim(IMGo,Gimg);
Sum = [S_MSE_f,S_PSNR_f,S_SSIM_f,S_MSE_o,S_PSNR_o,S_SSIM_o];
%------------------------------------------------------------------------------%
%subplot(1,2,1),imshow(Gimg),title('Original')
%subplot(1,2,2),imshow(IMGf),title('ART')
pseq = linspace(1,Gpix,Gpix);
hsec = Gpix/2;  %horizontal section
vsec = Gpix/2;  %vertical section
subplot(2,1,1),plot(pseq,Gimg(hsec,:),pseq,IMGf(hsec,:),'r');
subplot(2,1,2),plot(pseq,Gimg(:,vsec),pseq,IMGf(:,vsec),'r');
%imwrite(IMGf,strcat(FN,'.bmp'));
%imwrite(IMGo,strcat(FN,'NBG.bmp'));
%------------------------------------------------------------------------------%
toc  %结束计时
%------------------------------------------------------------------------------%