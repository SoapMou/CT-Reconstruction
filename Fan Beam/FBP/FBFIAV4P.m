%Fan-Beam FBP Interpolation A V4 Parfor
%射线:扇形束--探测器:等角
%投影:插值法--反投影:插值法--重建算法:FBP
%----------------------------------Description---------------------------------%
%图像点阵均匀分布在1~4象限,且间距为1;反投影时,将射线点阵平移到第一象限
%基于FBIAV3Par
%------------------------------------------------------------------------------%
tic  %开始计时
%----------------------------------Parameters----------------------------------%
Gsod  = 726;    %射线源到旋转中心的距离,单位为一个像素
Gfan  = 60;     %扇形束张角(角度制)
Ganum = 1080;   %旋转角度数量(射线源从3点钟方向逆时针旋转一周)
Gdnum = 768;    %探测器数量(探测器编号逆时针增加)
Gpix  = 512;    %图像尺寸(FOV的内接正方形,sod推荐为0.5*pix/(sind(fan/2)*sind(45))+2)
Gfnum = 1;      %滤波器编号(1=Ramp;2=Hanning;3=Hamming;4=Blackman;5=Bartlett;6=Triang;7=Bartlet-Hann)
Gflen = 16384;  %滤波器长度
%----------------------------------File Name-----------------------------------%
FN1 = strcat('FBFIAV4s',string(Gsod),'f',string(Gfan*10),'p',string(Gpix), ...
             'a',string(Ganum),'d',string(Gdnum),'f',string(Gfnum),'.bmp');
FN2 = strcat('FBFIAV4s',string(Gsod),'f',string(Gfan*10),'p',string(Gpix), ...
             'a',string(Ganum),'d',string(Gdnum),'f',string(Gfnum),'NBG.bmp');
%------------------------------------------------------------------------------%
Gimg  = phantom(Gpix);  %重建对象
Giba  = transpose(linspace(-0.5*Gfan,0.5*Gfan,Gdnum));  %扇形束各射线的初始角度(角度制)
GFBP  = zeros(Gpix*Gpix,1);  %FBP数据
maskB = zeros(Gpix*Gpix,1);  %maskB(背景区域≠0,目标区域=0)
%------------------------------Weighting Sequence------------------------------%
Gw = 1./cosd(Giba);  %加权序列(列向量)
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
%--------------------------------Image meshgrid--------------------------------%
Gms = (Gpix-1)/2;  %meshgrid shift(网格偏移量)
[GIX,GIY] = meshgrid(-Gms:1:Gms,-Gms:1:Gms);  %重建对象的XY坐标网格(上下颠倒)
%--------------------------Initial fan-beam meshgrid---------------------------%
Gfms = ceil(sqrt(2)*Gpix/2)+0.5;  %射线网格宽度/2
Gfw  = 2*Gfms+1;  %射线网格宽度(每条射线的采样点数量)
GFX  = repmat(-Gfms:1:Gfms,Gdnum,1);  %初始状态射线的X坐标网格
GFY  = tand(Giba).*(GFX-Gsod);  %初始状态射线的Y坐标网格
%-----------------------------FBP on each fan beam-----------------------------%
parfor a = 1:Ganum
    %--------将全局变量变为局部变量--------%
    sa   = (a-1)*360/Ganum;  %当前射线源的角度(角度制)
    pix  = Gpix;
    ms   = Gms;
    fw   = Gfw;
    FX   = GFX;
    FY   = GFY;
    %--------其他局部变量--------%
    mask = zeros(pix*pix,1);  %当前角度mask(背景区域=1,目标区域=0)
    FBP  = zeros(pix*pix,1);
    %--------Calculating projection--------%
    X = cosd(sa)*FX-sind(sa)*FY;  %当前角度射线X坐标网格(逆时针旋转)
    Y = sind(sa)*FX+cosd(sa)*FY;  %当前角度射线Y坐标网格(逆时针旋转)
    proj = interp2(GIX,GIY,Gimg,X,Y,'linear');  %使用二维插值法计算投影数据
    proj(isnan(proj)) = 0;  %将nan赋值0
    pdata = sum(proj,2);  %将proj每行数据累加得到每条射线的投影值
    %--------Weighting/Filtering/Backprojection--------% 
    wpdata  = pdata.*Gw;
    fftp    = fftshift(fft(wpdata,Gflen));  %fftshift使fft变为负高频-低频-高频形式,方便滤波
    ifftp   = real(ifft(ifftshift(fftp.*Gf)));  %滤波+反变换
    fwpdata = ifftp(1:Gdnum);
    for r = 1:Gdnum
        for c = 1:fw
            if abs(X(r,c)) <= ms && abs(Y(r,c)) <= ms
                xd = X(r,c)+ms-floor(X(r,c)+ms);  %当前插值点与左侧x坐标的差
                ct = ceil(X(r,c)+ms)+1;  %column top,当前插值点右侧列
                cb = floor(X(r,c)+ms)+1;  %column bottom,当前插值点左侧列
                yd = ceil(Y(r,c)+ms)-Y(r,c)-ms;  %当前插值点与下方y坐标的差(图像上下翻转后,下方y对应ceil)
                rt = floor(Y(r,c)+ms)+1;  %row top,当前插值点上方行
                rb = ceil(Y(r,c)+ms)+1;  %row bottom,当前插值点下方行
                FBP((ct-1)*pix+rt) = FBP((ct-1)*pix+rt)+xd*yd*fwpdata(r);  %给右上点赋值(注意:距离越近权重越大)
                FBP((cb-1)*pix+rt) = FBP((cb-1)*pix+rt)+(1-xd)*yd*fwpdata(r);  %给左上点赋值
                FBP((ct-1)*pix+rb) = FBP((ct-1)*pix+rb)+xd*(1-yd)*fwpdata(r);  %给右下点赋值
                FBP((cb-1)*pix+rb) = FBP((cb-1)*pix+rb)+(1-xd)*(1-yd)*fwpdata(r);  %给左下点赋值
                if pdata(r) == 0  %基于Bilinear方法确定背景区域
                    mask((ct-1)*pix+rt) = 1;  %右上点
                    mask((cb-1)*pix+rt) = 1;  %左上点
                    mask((ct-1)*pix+rb) = 1;  %右下点
                    mask((cb-1)*pix+rb) = 1;  %左下点
                end
            end
        end
    end
    maskB = maskB+mask;
    GFBP = GFBP+FBP;
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