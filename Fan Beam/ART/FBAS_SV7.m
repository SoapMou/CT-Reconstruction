%Fan-Beam ART Siddon Slow V7
%----------------------------------Description---------------------------------%
%系统矩阵:Siddon(即时计算系统矩阵,占用内存小,速度慢)
%重建算法:ART
%基础程序:FBAS_SV6
%------------------------------------------------------------------------------%
tic  %开始计时
%----------------------------------Parameters----------------------------------%
Gdet  = 'A';   %探测器排列方式(A=等角;L=等距)
Gpix  = 512;   %图像尺寸(单边像素数量)
Gsod  = 726;   %射线源-旋转中心距离(sod>0.5*pix/(sin(fan/2)*sin(45°)))
Gfan  = 60;    %扇形束张角(角度制)
Ganum = 1200;  %旋转角度数量(射线源沿逆时针方向旋转)
Gdnum = 801;   %探测器数量(探测器沿逆时针方向编号)
Gfidx = 1;     %滤波器编号(1=Ramp;2=Hanning;3=Hamming;4=Cosine)
Ginum = 6;     %迭代次数
RF    = 0.1;   %松弛因子(控制迭代的步长)
TV    = 0;     %TV正则(1=使用;0=不使用)
Mode  = 1;     %迭代重建模式(0=基于空图;1=基于FBP)
Gimg  = MAYOp512;  %重建对象
%-----------------------------------File Name----------------------------------%
FN = strcat('FBAS_SV7',Gdet,'p',string(Gpix),'s',string(Gsod),'f',string(Gfan*10),...
            'a',string(Ganum),'d',string(Gdnum),'f',string(Gfidx),'i',string(Ginum),...
            'r',string(RF),'m',string(Mode));
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
Gfilt = transpose(ifftshift(abs(linspace(-1,1,Gflen)).*win(Gfidx,:)));  %预先ifftshift
%-------------------------------FBP at each angle------------------------------%
parfor a = 1:Ganum
    %--------全局变量->局部变量--------%
    pix  = Gpix;
    dnum = Gdnum;
    img  = Gimg;
    sa   = (a-1)*360/Ganum;  %source angle(角度制)
    ba   = Giba+sa;  %beam angle(角度制)
    x0   = Gsod*cosd(sa);  %射线源的X轴坐标
    y0   = Gsod*sind(sa);  %射线源的Y轴坐标
    %--------其他局部变量--------%
    hp   = pix/2;
    FBP  = zeros(pix*pix,1);  %FBP数据
    mask = zeros(pix*pix,1);  %背景区域(背景区域=1,目标区域=0)
    SIDn = zeros(dnum,1);  %Siddon number(射线穿过像素的数量)
    SIDi = zeros(pix*2-1,dnum);  %Siddon index(射线穿过像素的编号)
    SIDw = ones(pix*2-1,dnum);  %Siddon weight(射线穿过像素的权重)(默认=1)
    PD   = zeros(dnum,1);  %投影数据
    COS  = abs(cosd(ba));
    %--------Projection--------%
    for d = 1:dnum
        k = tand(ba(d));  %射线的斜率
        if abs(k) == inf  %--------射线垂直于X轴--------%
            if abs(x0) < hp  %判断射线是否穿过图像
                col = ceil(x0)+hp;  %确定射线对应的列
                SIDn(d) = pix;  %记录像素的数量
                SIDi(1:pix,d) = (col-1)*pix+1:1:col*pix;  %记录像素的编号
                PD(d) = sum(img(:,col));  %记录射线的投影值
            end 
        elseif k == 0  %--------射线垂直于Y轴--------%
            if abs(y0) < hp  %判断射线是否穿过图像
                row = -floor(y0)+hp;  %确定射线对应的行
                SIDn(d) = pix;  %记录像素的数量
                SIDi(1:pix,d) = row:pix:(pix-1)*pix+row;  %记录像素的编号
                PD(d) = sum(img(row,:));  %记录射线的投影值
            end
        else  %--------射线处于其他角度--------%
            xT = (hp-y0)/k+x0;  %x top
            xB = (-hp-y0)/k+x0;  %x bottom
            if (xT<=-hp && xB<=-hp)+(xT>=hp && xB>=hp) ~= 0  %射线是否从图像左/右侧经过
                continue;
            end
            if k < 0  %射线与图像左/右交点X坐标
                X = [max(-hp,xT),min(hp,xB)];
            else
                X = [max(-hp,xB),min(hp,xT)];
            end
            Y = (X-x0)*k+y0;  %射线与图像左/右交点Y坐标
            xpos = unique([X,ceil(X(1)):floor(X(2)),((ceil(min(Y)):floor(max(Y)))-y0)/k+x0]);  %所有交点的X坐标
            SIDn(d) = length(xpos)-1;  %射线穿过像素的数量
            for n = 1:SIDn(d)
                row = hp-floor((xpos(n)+xpos(n+1))*k/2-x0*k+y0);  %当前像素的行
                col = hp+1+floor(xpos(n));  %当前像素的列(xpos从小到大,简化计算)
                SIDi(n,d) = (col-1)*pix+row;  %记录像素的编号
                SIDw(n,d) = (xpos(n+1)-xpos(n))/COS(d);  %记录像素的权重
                PD(d) = PD(d)+img(row,col)*SIDw(n,d);  %记录射线的投影值
            end
        end
    end
    %--------Filtering/Backprojection--------%
    fftPD = real(ifft(fft(PD,Gflen).*Gfilt));  %滤波
    fPD   = fftPD(1:dnum);
    for d = 1:dnum
        if SIDn(d) ~= 0  %判断射线是否穿过图像
            for n = 1:SIDn(d)  %反投影
                FBP(SIDi(n,d)) = FBP(SIDi(n,d))+fPD(d)*SIDw(n,d);
            end
            if PD(d) == 0  %记录背景区域
                mask(SIDi(1:SIDn(d),d)) = 1;
            end
        end
    end
    GFBP  = GFBP+FBP;  %将当前角度的FBP数据累加到GFBP中
    Gmask = Gmask+mask;  %将当前角度的mask累加到Gmask中
    sino(:,a) = PD;  %将当前角度的投影数据记录到Sinogram中
end
%---------------------------Iterative Reconstruction---------------------------%
FBPf = min(max(GFBP*pi/(2*Ganum),0),1);  %将FBP重建数据规格化并约束在[0,1]的范围
FBPo = FBPf.*(1-logical(Gmask));  %FBP重建数据(去除背景信息)
pix  = Gpix;
hp   = pix/2;
if Mode == 0
    ARTf = zeros(pix*pix,1);  %基于空图重建
elseif Mode == 1
    ARTf = FBPo;  %基于FBP结果重建
end
I = zeros(pix*2-1,1);  %index(射线穿过像素的编号)
W = zeros(pix*2-1,1);  %weight(射线穿过像素的权重)
for i = 1:Ginum
    for a = 1:Ganum
        sa  = (a-1)*360/Ganum;  %source angle(角度制)
        ba  = Giba+sa;  %beam angle(角度制)
        x0  = Gsod*cosd(sa);  %射线源的X轴坐标
        y0  = Gsod*sind(sa);  %射线源的Y轴坐标
        COS = abs(cosd(ba));
        %--------Update image at each beam--------%
        for d = 1:Gdnum  %每分析一条射线更新一次图像
            k = tand(ba(d));  %射线的斜率
            PDe = 0;  %投影值/投影值差值
            wt2 = 0;  %权重平方和
            if abs(k) == inf  %--------射线垂直于x轴--------%
                if abs(x0) < hp  %判断射线是否穿过图像
                    col = ceil(x0)+hp;  %确定射线对应的列
                    I(1:pix) = (col-1)*pix+1:1:col*pix;  %记录像素的编号
                    PDe = RF*(sino(d,a)-sum(ARTf(I(1:pix))))/pix;  %射线投影值差值
                    ARTf(I(1:pix)) = ARTf(I(1:pix))+PDe;  %更新图像
                end
            elseif k == 0  %--------射线垂直于y轴--------%
                if abs(y0) < hp  %判断射线是否穿过图像
                    row = -floor(y0)+hp;  %确定射线对应的行
                    I(1:pix) = row:pix:(pix-1)*pix+row;  %记录像素的编号
                    PDe = RF*(sino(d,a)-sum(ARTf(I(1:pix))))/pix;  %射线投影值差值
                    ARTf(I(1:pix)) = ARTf(I(1:pix))+PDe;  %更新图像
                end
            else  %--------射线处于其他角度--------%
                xT = (hp-y0)/k+x0;  %x top
                xB = (-hp-y0)/k+x0;  %x bottom
                if (xT<=-hp && xB<=-hp)+(xT>=hp && xB>=hp) ~= 0  %射线是否从图像左/右侧经过
                    continue;
                end
                if k < 0  %射线与图像左/右交点X坐标
                    X = [max(-hp,xT),min(hp,xB)];
                else
                    X = [max(-hp,xB),min(hp,xT)];
                end
                Y = (X-x0)*k+y0;  %射线与图像左/右交点Y坐标
                xpos = unique([X,ceil(X(1)):floor(X(2)),((ceil(min(Y)):floor(max(Y)))-y0)/k+x0]);  %所有交点的X坐标
                N = length(xpos)-1;  %记录像素的数量
                for n = 1:N
                    row = hp-floor((xpos(n)+xpos(n+1))*k/2-x0*k+y0);  %当前像素的行
                    col = hp+1+floor(xpos(n));  %当前像素的列(xpos从小到大,简化计算)
                    I(n) = (col-1)*pix+row;  %记录像素的编号
                    W(n) = (xpos(n+1)-xpos(n))/COS(d);  %记录像素的权重
                    PDe = PDe+ARTf((col-1)*pix+row)*W(n);  %记录射线的投影值
                    wt2 = wt2+W(n)^2;
                end
                PDe = RF*(sino(d,a)-PDe)/wt2;  %射线投影值差值
                for n = 1:N  %更新图像
                    ARTf(I(n)) = ARTf(I(n))+W(n)*PDe;
                end
            end
        end
    end
    ARTf = min(max(ARTf,0),1);  %约束<0和>1的值
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