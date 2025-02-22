%Fan-Beam FBP Interpolation V4 Parfor
%----------------------------------Description---------------------------------%
%系统矩阵:Interpolation
%重建算法:FBP
%基础程序:FBIAV3Par
%------------------------------------------------------------------------------%
tic  %开始计时
%----------------------------------Parameters----------------------------------%
Gdet  = 'A';   %探测器排列方式(A=等角;L=等距)
Gpix  = 512;   %图像尺寸(单边像素数量)
Gsod  = 726;   %射线源-旋转中心距离(sod>0.5*pix/(sin(fan/2)*sin(45°)))
Gfan  = 60;    %扇形束张角(角度制)
Ganum = 1152;  %旋转角度数量(射线源沿逆时针方向旋转)
Gdnum = 769;   %探测器数量(探测器沿逆时针方向编号)
Gfidx = 1;     %滤波器编号(1=Ramp;2=Hanning;3=Hamming;4=Cosine)
Gimg  = phantom(Gpix);  %重建对象
%----------------------------------File Name-----------------------------------%
FN = strcat('FBFI_V4',Gdet,'p',string(Gpix),'s',string(Gsod),'f',string(Gfan*10),...
            'a',string(Ganum),'d',string(Gdnum),'f',string(Gfidx));
%------------------------------------------------------------------------------%
GFBP  = zeros(Gpix*Gpix,1);  %FBP数据
Gmask = zeros(Gpix*Gpix,1);  %背景区域(背景区域>0,目标区域=0)
if strcmp(Gdet,'A')  %射线的初始角度(等角排布)
    Giba = transpose(linspace(-Gfan/2,Gfan/2,Gdnum));  
elseif strcmp(Gdet,'L')  %射线的初始角度(等距排布)
    Giba = transpose(atand(linspace(tand(-Gfan/2),tand(Gfan/2),Gdnum)));
end
%------------------------------------Filter------------------------------------%
Gflen = 4096;  %滤波器长度
win   = ones(10,Gflen);  %窗口函数
win(2,:) = hann(Gflen);
win(3,:) = hamming(Gflen);
win(4,:) = cos(pi*linspace(-1,1,Gflen)/2);
Gfilt = transpose(abs(linspace(-1,1,Gflen)).*win(Gfidx,:));
%--------------------------------Image meshgrid--------------------------------%
Gms = (Gpix-1)/2;  %meshgrid shift(图像网格偏移量)
[GIX,GIY] = meshgrid(-Gms:1:Gms,Gms:-1:-Gms);  %重建对象的XY坐标网格
%---------------------------Initial fan-beam meshgrid--------------------------%
Gfms = ceil(sqrt(2)*Gpix/2)+0.5;  %射线网格偏移量
Gfw  = 2*Gfms+1;  %射线网格宽度(每条射线的采样点数量)
GFX  = (repmat(-Gfms:1:Gfms,Gdnum,1)-Gsod).*cosd(Giba)+Gsod;  %初始状态射线的X坐标网格
GFY  = tand(Giba).*(GFX-Gsod);  %初始状态射线的Y坐标网格
%-------------------------------FBP at each angle------------------------------%
parfor a = 1:Ganum
    %--------全局变量->局部变量--------%
    pix  = Gpix;
    ms   = Gms;
    fw   = Gfw;
    FX   = GFX;
    FY   = GFY;
    sa   = (a-1)*360/Ganum;  %source angle(角度制)
    %--------其他局部变量--------%
    FBP  = zeros(pix*pix,1);  %FBP数据
    mask = zeros(pix*pix,1);  %背景区域(背景区域=1,目标区域=0)
    %--------Projection of each beam--------%
    X = cosd(sa)*FX-sind(sa)*FY;  %射线的X坐标网格(逆时针旋转)
    Y = sind(sa)*FX+cosd(sa)*FY;  %射线的Y坐标网格(逆时针旋转)
    proj = interp2(GIX,GIY,Gimg,X,Y,'linear');  %使用线性插值法计算投影数据
    proj(isnan(proj)) = 0;  %将nan赋值0
    PD = sum(proj,2);  %将proj每行数据累加得到每条射线的投影值
    %--------Filtering/Backprojection--------%
    fftp  = fftshift(fft(PD,Gflen));  %fftshift将fft变为负高频-低频-高频形式
    ifftp = real(ifft(ifftshift(fftp.*Gfilt)));  %滤波+傅里叶反变换
    fPD   = ifftp(1:Gdnum);
    for r = 1:Gdnum
        for c = 1:fw
            if abs(X(r,c)) <= ms && abs(Y(r,c)) <= ms
                dx = X(r,c)+ms-floor(X(r,c)+ms);  %插值点与左侧X坐标的差
                ct = ceil(X(r,c)+ms)+1;  %column top,插值点右侧对应图像的列
                cb = floor(X(r,c)+ms)+1;  %column bottom,插值点左侧对应图像的列
                dy = Y(r,c)+ms-floor(Y(r,c)+ms);  %插值点与下方Y坐标的差
                rt = pix-ceil(Y(r,c)+ms);  %row top,插值点上方对应图像的行
                rb = pix-floor(Y(r,c)+ms);  %row bottom,插值点下方对应图像的行
                FBP((ct-1)*pix+rt) = FBP((ct-1)*pix+rt)+dx*dy*fPD(r);  %右上
                FBP((cb-1)*pix+rt) = FBP((cb-1)*pix+rt)+(1-dx)*dy*fPD(r);  %左上
                FBP((ct-1)*pix+rb) = FBP((ct-1)*pix+rb)+dx*(1-dy)*fPD(r);  %右下
                FBP((cb-1)*pix+rb) = FBP((cb-1)*pix+rb)+(1-dx)*(1-dy)*fPD(r);  %左下
                if PD(r) == 0  %投影值=0,记录背景区域
                    mask((ct-1)*pix+rt) = 1;  %右上点
                    mask((cb-1)*pix+rt) = 1;  %左上点
                    mask((ct-1)*pix+rb) = 1;  %右下点
                    mask((cb-1)*pix+rb) = 1;  %左下点
                end
            end
        end
    end
    GFBP = GFBP+FBP;  %将当前角度的FBP数据累加到GFBP中
    Gmask = Gmask+mask;  %将当前角度的mask累加到Gmask中
end
%-----------------------Normalization and Reconstruction-----------------------%
FBPf = min(max(GFBP*pi/(2*Ganum),0),1);  %将FBP重建数据规格化并约束在[0,1]的范围
FBPo = FBPf.*(1-logical(Gmask));  %FBP重建数据(去除背景信息)
IMGf = reshape(FBPf,Gpix,Gpix);  %FBP重建图像(完整)
IMGo = reshape(FBPo,Gpix,Gpix);  %FBP重建图像(去除背景信息)
%------------------------Evaluate Reconstruction Quality-----------------------%
S_MSE_f  = mean((FBPf-Gimg(:)).^2);
S_MSE_o  = mean((FBPo-Gimg(:)).^2);
S_PSNR_f = psnr(IMGf,Gimg);
S_PSNR_o = psnr(IMGo,Gimg);
S_SSIM_f = ssim(IMGf,Gimg);
S_SSIM_o = ssim(IMGo,Gimg);
Sum = [S_MSE_f,S_PSNR_f,S_SSIM_f,S_MSE_o,S_PSNR_o,S_SSIM_o];
%------------------------------------------------------------------------------%
%subplot(1,2,1),imshow(Gimg),title('Original')
%subplot(1,2,2),imshow(IMGf),title('FBP')
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