%Fan-Beam FBP Siddon A V9 Parfor
%射线:扇形束--探测器:等角
%投影:Siddon--反投影:Siddon--重建算法:FBP
%----------------------------------Description---------------------------------%
%图像均匀分布在1~4象限
%基于FBV8RPar
%------------------------------------------------------------------------------%
tic  %开始计时
%----------------------------------Parameters----------------------------------%
Gsod  = 726;    %射线源到旋转中心的距离,单位为一个像素
Gfan  = 60;     %扇形束张角(角度制)
Ganum = 1200;   %旋转角度数量(射线源从3点钟方向逆时针旋转一周)
Gdnum = 801;    %探测器数量(探测器编号逆时针增加)
Gpix  = 512;    %图像尺寸(FOV的内接正方形,sod推荐为0.5*pix/(sind(fan/2)*sind(45))+2)
Gfnum = 1;      %滤波器编号(1=Ramp;2=Hanning;3=Hamming;4=Blackman;5=Bartlett;6=Triang;7=Bartlet-Hann)
Gflen = 16384;  %滤波器长度
%----------------------------------File Name-----------------------------------%
FN1 = strcat('FBFSAV9s',string(Gsod),'f',string(Gfan*10),'p',string(Gpix), ...
             'a',string(Ganum),'d',string(Gdnum),'f',string(Gfnum),'.bmp');
FN2 = strcat('FBFSAV9s',string(Gsod),'f',string(Gfan*10),'p',string(Gpix), ...
             'a',string(Ganum),'d',string(Gdnum),'f',string(Gfnum),'NBG.bmp');
%------------------------------------------------------------------------------%
Gimg  = MAYOp512;  %重建对象
Giseq = Gimg(:);  %重建对象的列向量形式
Giba  = transpose(linspace(-0.5*Gfan,0.5*Gfan,Gdnum));  %扇形束各射线的初始角度(角度制)
GFBP  = zeros(Gpix*Gpix,1);  %FBP数据
maskB = zeros(Gpix*Gpix,1);  %maskB(背景区域≠0,目标区域=0)
%------------------------------Weighting Sequence------------------------------%
Gw = 1./cosd(Giba);  %加权序列
%------------------------------------Filter------------------------------------%
filter = zeros(Gflen,10);
ramp = transpose(abs(linspace(-1,1,Gflen)));  %斜坡滤波器
filter(:,1) = ramp;                      %1:Ramp
filter(:,2) = ramp.*hann(Gflen);         %2:Hanning
filter(:,3) = ramp.*hamming(Gflen);      %3:Hamming
filter(:,4) = ramp.*blackman(Gflen);     %4:Blackman
filter(:,5) = ramp.*bartlett(Gflen);     %5:Bartlett
filter(:,6) = ramp.*triang(Gflen);       %6:Triang
filter(:,7) = ramp.*barthannwin(Gflen);  %7:Bartlet-Hann
Gf = filter(:,Gfnum);  %当前选择的滤波器
%-----------------------------FBP on each fan beam-----------------------------%
parfor a = 1:Ganum
    %--------将全局变量变为局部变量--------%
    dnum  = Gdnum;
    pix   = Gpix;
    img   = Gimg;
    iseq  = Giseq;
    sa    = (a-1)*360/Ganum;  %当前射线源的角度(角度制)
    ba    = Giba+sa;  %当前扇形束各射线的角度(角度制)
    %--------其他局部变量--------%
    hpix  = 0.5*pix;
    x0    = Gsod*cosd(sa);  %当前角度射线源的横坐标
    y0    = Gsod*sind(sa);  %当前角度射线源的纵坐标
    FBP   = zeros(pix*pix,1);  %当前角度FBP数据
    pdata = zeros(dnum,1);  %当前角度投影数据
    cmat  = zeros(pix*2-1,dnum);  %坐标矩阵(每列记录每条射线穿过像素的位置)
    bmat  = zeros(pix*2-1,dnum);  %权重矩阵(每列记录每条射线穿过像素的权重)
    nseq  = zeros(dnum,1);  %计数向量(记录每条射线穿过像素的数量)
    mask  = zeros(pix*pix,1);  %当前角度mask(背景区域=1,目标区域=0)
    %--------Calculate Coordinate/Bevel/Projection--------%
    for d = 1:dnum
        k = tand(ba(d));  %当前射线的斜率
        %--------当射线垂直于X轴时--------%
        if abs(k) == inf
            if abs(x0) <= hpix  %判断射线是否穿过图像
                col = ceil(x0)+hpix;  %确定射线对应的列
                if col == 0
                    col = col+1;
                end
                pdata(d) = sum(img(:,col));  %计算投影数据
                cmat(1:pix,d) = (col-1)*pix+1:1:col*pix;  %记录射线穿过像素的位置
                bmat(1:pix,d) = ones(pix,1);  %记录射线穿过像素的权重
                nseq(d) = pix;  %射线穿过了pix个像素
                if pdata(d) == 0  %投影值=0时,将射线的路径记录到mask中
                    mask(cmat(1:pix,d)) = 1;
                end
            end
        %--------当射线垂直于Y轴时--------%
        elseif k == 0
            if abs(y0) <= hpix  %判断射线是否穿过图像
                row = -floor(y0)+hpix;  %确定射线对应的行
                if row == 0
                    row = row+1;
                end
                pdata(d) = sum(img(row,:));  %计算投影数据
                cmat(1:pix,d) = row:pix:(pix-1)*pix+row;  %记录射线穿过像素的位置
                bmat(1:pix,d) = ones(pix,1);  %记录射线穿过像素的权重
                nseq(d) = pix;  %射线穿过了pix个像素
                if pdata(d) == 0  %投影值=0时,将射线的路径记录到mask中
                    mask(cmat(1:pix,d)) = 1;
                end
            end
        %--------当射线处于其他角度时--------%
        else
            partx = x0-y0/k;  %x=y/k+partx(已知y,求x)
            party = y0-x0*k;  %y=x*k+party(已知x,求y)
            xtop = hpix/k+partx;
            xbutt = -hpix/k+partx;
            Rsec = xtop>hpix && xbutt>hpix;  %判断射线是否从图像右侧经过
            Lsec = xtop<-hpix && xbutt<-hpix;  %判断射线是否从图像左侧经过
            if Rsec+Lsec ~= 0  %射线与图像无交点时跳过本次循环
                continue;
            end
            cos = abs(cosd(ba(d)));
            %--------k<0时,射线方向朝右下--------%
            if k < 0
                xL = max(-hpix,xtop);  %计算射线与图像左交点X坐标
                xR = min(hpix,xbutt);  %计算射线与图像右交点X坐标
            %--------k>0时,射线方向朝右上--------%
            else
                xL = max(-hpix,xbutt);  %计算射线与图像左交点X坐标
                xR = min(hpix,xtop);  %计算射线与图像右交点X坐标
            end
            yLR = sort([xL*k+party,xR*k+party]);  %计算左右交点Y坐标并升序排列
            %--------计算投影值和投影矩阵--------%
            seqy = floor(yLR(1))+1:1:ceil(yLR(2))-1;  %中间所有整数Y坐标(不包括左右交点)
            allx = unique([xL,floor(xL)+1:1:ceil(xR)-1,xR,seqy/k+partx]);  %射线与线框相交的所有X坐标(使用unique去掉重复项并排序)
            ally = allx*k+party;  %计算射线与线框相交的所有Y坐标
            xlen = length(allx);
            nseq(d) = xlen-1;  %射线穿过了xlen-1个像素
            rnow = hpix-floor((ally(1:xlen-1)+ally(2:xlen))*0.5);  %射线穿过像素的行数
            cnow = hpix+1+floor(allx(1:xlen-1));  %射线穿过像素的列数(x固定为从小到大,简化计算)
            cord = (cnow-1)*pix+rnow;  %射线穿过像素的位置
            bevel = transpose((allx(2:xlen)-allx(1:xlen-1))/cos);  %射线穿过像素的权重
            cmat(1:xlen-1,d) = cord;
            bmat(1:xlen-1,d) = bevel;
            pdata(d) = sum(iseq(cord).*bevel);  %记录当前射线的投影值
            if pdata(d) == 0  %投影值=0时,将射线的路径记录到mask中
                mask(cord) = 1;
            end
        end
    end
    %--------Weighting/Filtering/Backprojection--------%
    wpdata  = pdata;  %实验表明Siddon算法不加权的重建效果更好
    fftp    = fftshift(fft(wpdata,Gflen));  %fftshift使fft变为负高频-低频-高频形式,方便滤波
    ifftp   = real(ifft(ifftshift(fftp.*Gf)));  %滤波+反变换
    fwpdata = ifftp(1:dnum);
    for d = 1:dnum
        cord = cmat(1:nseq(d),d);  %当前射线穿过像素的位置
        FBP(cord) = FBP(cord)+fwpdata(d)*bmat(1:nseq(d),d);
    end
    GFBP = GFBP+FBP;  %将当前角度的FBP数据累加到GFBP中
    maskB = maskB+mask;  %将当前角度的mask累加到maskB中
end
%-------------------------------Process FBP Data-------------------------------%
GFBP = max(GFBP,0);  %将<0的值归零
maskO = 1-logical(maskB);  %maskO(目标区域=1;背景区域=0)
oseq = GFBP(maskO==1);  %记录目标区域的FBP数据
omax = max(oseq);  %目标区域最大值
omin = min(oseq);  %目标区域最小值
FBPn = (GFBP-omin)/(omax-omin);  %以目标区域为准,对整个图像进行规格化
FBPn = max(FBPn,0);  %再次将<0的值归零
FBPnbg = FBPn.*maskO;  %去除背景区域的规格化数据
%------------------------------Reconstruct Image-------------------------------%
FBPimg = reshape(FBPn,Gpix,Gpix);  %还原普通重建图像
FBPnbgimg = reshape(FBPnbg,Gpix,Gpix);  %还原去除背景区域的图像
%------------------------Evaluate Reconstruction Quality-----------------------%
S_MSE1  = mean((FBPn-Gimg(:)).^2);
S_MSE2  = mean((FBPnbg-Gimg(:)).^2);
S_PSNR1 = psnr(FBPimg,Gimg);
S_PSNR2 = psnr(FBPnbgimg,Gimg);
S_SSIM1 = ssim(FBPimg,Gimg);
S_SSIM2 = ssim(FBPnbgimg,Gimg);
Sum = [S_MSE1,S_PSNR1,S_SSIM1,S_MSE2,S_PSNR2,S_SSIM2];
%-----------------------------Plot and Export Image----------------------------%
%subplot(1,2,1), imshow(Gimg), title('Original')
%subplot(1,2,2), imshow(FBPimg), title('FBP')
x = linspace(1,Gpix,Gpix);
hsec = Gpix/2;
vsec = Gpix/2;
subplot(2,1,1), plot(x,Gimg(hsec,:),x,FBPimg(hsec,:),'r');
subplot(2,1,2), plot(x,Gimg(:,vsec),x,FBPimg(:,vsec),'r');
%imwrite(FBPimg,FN1);
%imwrite(FBPnbgimg,FN2);
%------------------------------------------------------------------------------%
toc  %停止计时
%------------------------------------------------------------------------------%