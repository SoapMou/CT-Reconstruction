%Fan-Beam ART Siddon A V7
%射线:扇形束--探测器:等角
%投影:Siddon--反投影:Siddon--重建算法:ART
%----------------------------------Description---------------------------------%
%ART:对比射线的当前投影值和理论投影值并将差值反投影;差值需乘以"1/bevel平方和"的权重和relax
%基于FBASAV6;不使用parfor;预先计算完整系统矩阵(占用内存大,速度快)
%------------------------------------------------------------------------------%
tic  %开始计时
%----------------------------------Parameters----------------------------------%
sod   = 726;    %射线源到旋转中心的距离,单位为一个像素
fan   = 60;     %扇形束张角(角度制)
anum  = 1200;   %旋转角度数量(射线源从3点钟方向逆时针旋转一周)
dnum  = 801;    %探测器数量(探测器编号逆时针增加)
pix   = 512;    %图像尺寸(FOV的内接正方形,sod推荐为0.5*pix/(sind(fan/2)*sind(45))+2)
fnum  = 1;      %滤波器编号(1=Ramp;2=Hanning;3=Hamming;4=Blackman;5=Bartlett;6=Triang;7=Bartlet-Hann)
flen  = 16384;  %滤波器长度
inum  = 6;      %迭代次数
relax = 0.1;    %松弛因子(每次迭代的步长)
TV    = 0;      %是否使用TV正则(1=使用;0=不使用)
Mode  = 1;      %0=基于空图重建;1=基于FBP结果重建
Range = 1;      %0=每轮迭代后约束<0的值;1=每轮迭代后约束<0和>1的值
%----------------------------------File Name-----------------------------------%
FN1 = strcat('FBASAV7s',string(sod),'f',string(fan*10),'p',string(pix),'a',string(anum), ...
             'd',string(dnum),'f',string(fnum),'i',string(inum),'r',string(relax),'.bmp');
FN2 = strcat('FBASAV7s',string(sod),'f',string(fan*10),'p',string(pix),'a',string(anum), ...
             'd',string(dnum),'f',string(fnum),'i',string(inum),'r',string(relax),'NBG.bmp');
%------------------------------------------------------------------------------%
img   = phantom(pix);  %重建对象
iseq  = img(:);  %重建对象的列向量形式
iba   = transpose(linspace(-0.5*fan,0.5*fan,dnum));  %扇形束各射线的初始角度(角度制)
FBP   = zeros(pix*pix,1);  %FBP数据
maskB = zeros(pix*pix,1);  %maskB(背景区域≠0,目标区域=0)
sino  = zeros(dnum,anum);  %sinogram(投影数据)
%---------------------------------System Matrix--------------------------------%
cmat  = zeros(pix*2-1,dnum*anum);  %坐标矩阵(每列记录每条射线穿过像素的位置)
bmat  = zeros(pix*2-1,dnum*anum);  %权重矩阵(每列记录每条射线穿过像素的权重)
nmat  = zeros(dnum,anum);  %计数矩阵(记录每条射线穿过像素的数量)
%------------------------------Weighting Sequence------------------------------%
w = 1./cosd(iba);  %加权序列
%------------------------------------Filter------------------------------------%
filter = zeros(flen,10);
ramp = transpose(abs(linspace(-1,1,flen)));  %斜坡滤波器
filter(:,1) = ramp;                     %1:Ramp
filter(:,2) = ramp.*hann(flen);         %2:Hanning
filter(:,3) = ramp.*hamming(flen);      %3:Hamming
filter(:,4) = ramp.*blackman(flen);     %4:Blackman
filter(:,5) = ramp.*bartlett(flen);     %5:Bartlett
filter(:,6) = ramp.*triang(flen);       %6:Triang
filter(:,7) = ramp.*barthannwin(flen);  %7:Bartlet-Hann
f = filter(:,fnum);  %当前选择的滤波器
%----------------------------Calculate System Matrix---------------------------%
for a = 1:anum
    sa    = (a-1)*360/anum;  %当前射线源的角度(角度制)
    ba    = iba+sa;  %当前扇形束各射线的角度(角度制)
    hpix  = 0.5*pix;
    x0    = sod*cosd(sa);  %当前角度射线源的横坐标
    y0    = sod*sind(sa);  %当前角度射线源的纵坐标
    for d = 1:dnum
        k = tand(ba(d));  %当前射线的斜率
        %--------当射线垂直于X轴时--------%
        if abs(k) == inf
            if abs(x0) <= hpix  %判断射线是否穿过图像
                col = ceil(x0)+hpix;  %确定射线对应的列
                cmat(1:pix,(a-1)*dnum+d) = (col-1)*pix+1:1:col*pix;  %记录射线穿过像素的位置
                bmat(1:pix,(a-1)*dnum+d) = ones(pix,1);  %记录射线穿过像素的权重
                nmat(d,a) = pix;  %射线穿过了pix个像素
            end
        %--------当射线垂直于Y轴时--------%
        elseif k == 0
            if abs(y0) <= hpix  %判断射线是否穿过图像
                row = -floor(y0)+hpix;  %确定射线对应的行
                cmat(1:pix,(a-1)*dnum+d) = row:pix:(pix-1)*pix+row;  %记录射线穿过像素的位置
                bmat(1:pix,(a-1)*dnum+d) = ones(pix,1);  %记录射线穿过像素的权重
                nmat(d,a) = pix;  %射线穿过了pix个像素
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
            rnow = hpix-floor((ally(1:xlen-1)+ally(2:xlen))*0.5);  %射线穿过像素的行数
            cnow = hpix+1+floor(allx(1:xlen-1));  %射线穿过像素的列数(x固定为从小到大,简化计算)
            cmat(1:xlen-1,(a-1)*dnum+d) = (cnow-1)*pix+rnow;  %射线穿过像素的位置
            bmat(1:xlen-1,(a-1)*dnum+d) = transpose((allx(2:xlen)-allx(1:xlen-1))/cos);  %射线穿过像素的权重
            nmat(d,a) = xlen-1;  %射线穿过了xlen-1个像素
        end
    end
end
%-----------------------------FBP on each fan beam-----------------------------%
for a = 1:anum
    pdata = zeros(dnum,1);  %当前角度投影数据
    %--------Projection--------%
    for d = 1:dnum
        cord = cmat(1:nmat(d,a),(a-1)*dnum+d);  %当前射线穿过像素的位置
        pdata(d) = sum(iseq(cord).*bmat(1:nmat(d,a),(a-1)*dnum+d));  %记录当前射线的投影值
        if pdata(d) == 0  %投影值=0时,将射线的路径记录到mask中
            maskB(cord) = 1;
        end
    end
    sino(:,a) = pdata;  %将当前角度的投影数据记录到sinogram中
    %--------Weighting/Filtering/Backprojection--------%
    wpdata  = pdata;  %实验表明Siddon算法不加权的重建效果更好
    fftp    = fftshift(fft(wpdata,flen));  %fftshift使fft变为负高频-低频-高频形式,方便滤波
    ifftp   = real(ifft(ifftshift(fftp.*f)));  %滤波+反变换
    fwpdata = ifftp(1:dnum);
    for d = 1:dnum
        cord = cmat(1:nmat(d,a),(a-1)*dnum+d);  %当前射线穿过像素的位置
        FBP(cord) = FBP(cord)+fwpdata(d)*bmat(1:nmat(d,a),(a-1)*dnum+d);
    end
end
%-------------------------------Process FBP Data-------------------------------%
FBP = max(FBP,0);  %将<0的值归零
maskO = 1-maskB;  %maskO(目标区域=1;背景区域=0)
oseq = FBP(maskO==1);  %记录目标区域的FBP数据
omax = max(oseq);  %目标区域最大值
omin = min(oseq);  %目标区域最小值
FBPn = (FBP-omin)/(omax-omin);  %以目标区域为准,对整个图像进行规格化
FBPn = max(FBPn,0);  %再次将<0的值归零
%---------------------------Iterative Reconstruction---------------------------%
if Mode == 0
    cimg = zeros(pix*pix,1);  %基于空图重建
else
    cimg = FBPn;  %基于FBP结果重建
end
for i = 1:inum
    for a = 1:anum
        for d = 1:dnum
            cord = cmat(1:nmat(d,a),(a-1)*dnum+d);  %当前射线穿过像素的位置
            if sino(d,a) == 0  %当原始投影值=0时,当前路径的像素赋值0
                cimg(cord) = 0;
            else
                bevel = bmat(1:nmat(d,a),(a-1)*dnum+d);  %当前射线穿过像素的权重
                proj = sum(cimg(cord).*bevel);  %当前射线的投影值
                cimg(cord) = cimg(cord)+relax*bevel*(sino(d,a)-proj)/sum(bevel.^2);  %更新img
            end
        end
    end
    if Range == 0
        cimg = max(cimg,0);  %每轮迭代后约束<0的值
    else
        cimg = min(max(cimg,0),1);  %每轮迭代后约束<0和>1的值
    end
    %--------TV Regularization--------%
    if TV == 1
        img2 = reshape(cimg,pix,pix);
        [gradx,grady] = gradient(img2);
        gmag = sqrt(gradx.^2+grady.^2+1e-8);  %gradient magnitude
        div = divergence(gradx./gmag,grady./gmag);
        img2 = img2+0.0025*div;
        cimg = img2(:);
    end
end
%------------------------------Reconstruct Image-------------------------------%
FBPnbg = cimg.*maskO;  %获取去除背景区域的数据
FBPimg = reshape(cimg,pix,pix);  %还原迭代重建图像
FBPnbgimg = reshape(FBPnbg,pix,pix);  %还原去除背景区域的迭代重建图像
%-----------------------Evaluate Reconstruction Quality------------------------%
S_MSE1  = mean((cimg-iseq).^2);
S_MSE2  = mean((FBPnbg-iseq).^2);
S_PSNR1 = psnr(FBPimg,img);
S_PSNR2 = psnr(FBPnbgimg,img);
S_SSIM1 = ssim(FBPimg,img);
S_SSIM2 = ssim(FBPnbgimg,img);
Sum = [S_MSE1,S_PSNR1,S_SSIM1,S_MSE2,S_PSNR2,S_SSIM2];
%----------------------------Plot and Export Images----------------------------%
subplot(1,2,1), imshow(img), title('Original')
subplot(1,2,2), imshow(FBPimg), title('ART')
x = linspace(1,pix,pix);
hsec = pix/2;
vsec = pix/2;
%subplot(2,1,1), plot(x,img(hsec,:),x,FBPimg(hsec,:),'r');
%subplot(2,1,2), plot(x,img(:,vsec),x,FBPimg(:,vsec),'r');
%imwrite(FBPimg,FN1);
%imwrite(FBPnbgimg,FN2);
%------------------------------------------------------------------------------%
toc  %停止计时
%------------------------------------------------------------------------------%