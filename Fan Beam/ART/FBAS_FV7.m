%Fan-Beam ART Siddon Fast V7
%----------------------------------Description---------------------------------%
%系统矩阵:Siddon(预先计算系统矩阵,占用内存大,速度快)
%重建算法:ART
%基础程序:FBAS_FV6
%------------------------------------------------------------------------------%
tic  %开始计时
%----------------------------------Parameters----------------------------------%
det   = 'A';   %探测器排列方式(A=等角;L=等距)
pix   = 512;   %图像尺寸(单边像素数量)
sod   = 726;   %射线源-旋转中心距离(sod>0.5*pix/(sin(fan/2)*sin(45°)))
fan   = 60;    %扇形束张角(角度制)
anum  = 1200;  %旋转角度数量(射线源沿逆时针方向旋转)
dnum  = 801;   %探测器数量(探测器沿逆时针方向编号)
fidx  = 1;     %滤波器编号(1=Ramp;2=Hanning;3=Hamming;4=Cosine)
inum  = 8;     %迭代次数
RF    = 0.1;   %松弛因子(控制迭代的步长)
TV    = 0;     %TV正则(1=使用;0=不使用)
Mode  = 1;     %迭代重建模式(0=基于空图;1=基于FBP)
img   = phantom(pix);  %重建对象
%----------------------------------File Name-----------------------------------%
FN = strcat('FBAS_FV7',det,'p',string(pix),'s',string(sod),'f',string(fan*10),...
            'a',string(anum),'d',string(dnum),'f',string(fidx),'i',string(inum),...
            'r',string(RF),'m',string(Mode));
%------------------------------------------------------------------------------%
iseq  = img(:);  %重建对象(列向量)
hp    = pix/2;
FBP   = zeros(pix*pix,1);  %FBP数据
mask  = zeros(pix*pix,1);  %背景区域(背景区域>0,目标区域=0)
if strcmp(det,'A')  %射线的初始角度(等角排布)
    iba = transpose(linspace(-fan/2,fan/2,dnum));  
elseif strcmp(det,'L')  %射线的初始角度(等距排布)
    iba = transpose(atand(linspace(tand(-fan/2),tand(fan/2),dnum)));
end
SIDn  = zeros(dnum,anum);  %Siddon number(射线穿过像素的数量)
SIDi  = zeros(pix*2-1,dnum*anum);  %Siddon index(射线穿过像素的编号)
SIDw  = ones(pix*2-1,dnum*anum);  %Siddon weight(射线穿过像素的权重)(默认=1)
sino  = zeros(dnum,anum);  %Sinogram
%------------------------------------Filter------------------------------------%
flen  = 4096;  %滤波器长度
win   = ones(10,flen);  %窗口函数
win(2,:) = hann(flen);
win(3,:) = hamming(flen);
win(4,:) = cos(pi*linspace(-1,1,flen)/2);
filt = transpose(ifftshift(abs(linspace(-1,1,flen)).*win(fidx,:)));  %预先ifftshift
%-------------------------------FBP at each angle------------------------------%
for a = 1:anum
    sa = (a-1)*360/anum;  %source angle(角度制)
    ba = iba+sa;  %beam angle(角度制)
    x0 = sod*cosd(sa);  %射线源的X轴坐标
    y0 = sod*sind(sa);  %射线源的Y轴坐标
    COS = abs(cosd(ba));
    for d = 1:dnum
        k = tand(ba(d));  %射线的斜率
        if abs(k) == inf  %--------射线垂直于X轴--------%
            if abs(x0) < hp  %判断射线是否穿过图像
                col = ceil(x0)+hp;  %确定射线对应的列
                SIDn(d,a) = pix;  %记录像素的数量
                SIDi(1:pix,(a-1)*dnum+d) = (col-1)*pix+1:1:col*pix;  %记录像素的编号
                sino(d,a) = sum(img(:,col));  %记录射线的投影值
            end
        elseif k == 0  %--------射线垂直于Y轴--------%
            if abs(y0) < hp  %判断射线是否穿过图像
                row = -floor(y0)+hp;  %确定射线对应的行
                SIDn(d,a) = pix;  %记录像素的数量
                SIDi(1:pix,(a-1)*dnum+d) = row:pix:(pix-1)*pix+row;  %记录像素的编号
                sino(d,a) = sum(img(row,:));  %记录射线的投影值
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
            SIDn(d,a) = length(xpos)-1;  %记录像素的数量
            for n = 1:SIDn(d,a)
                row = hp-floor((xpos(n)+xpos(n+1))*k/2-x0*k+y0);  %当前像素的行
                col = hp+1+floor(xpos(n));  %当前像素的列(xpos从小到大,简化计算)
                SIDi(n,(a-1)*dnum+d) = (col-1)*pix+row;  %记录像素的编号
                SIDw(n,(a-1)*dnum+d) = (xpos(n+1)-xpos(n))/COS(d);  %记录像素的权重
                sino(d,a) = sino(d,a)+iseq((col-1)*pix+row)*SIDw(n,(a-1)*dnum+d);  %记录射线的投影值
            end
        end
    end
    %--------Filtering/Backprojection--------%
    fftPD = real(ifft(fft(sino(:,a),flen).*filt));  %滤波
    fPD   = fftPD(1:dnum);
    for d = 1:dnum
        if SIDn(d,a) ~= 0  %判断射线是否穿过图像
            for n = 1:SIDn(d,a)  %反投影
                FBP(SIDi(n,(a-1)*dnum+d)) = FBP(SIDi(n,(a-1)*dnum+d))+fPD(d)*SIDw(n,(a-1)*dnum+d);
            end
            if sino(d,a) == 0  %记录背景区域
                mask(SIDi(1:SIDn(d,a),(a-1)*dnum+d)) = 1;
            end
        end
    end
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
            if SIDn(d,a) ~= 0  %判断射线是否穿过图像
                cd = (a-1)*dnum+d;
                if sino(d,a) == 0  %原始投影值=0,射线穿过像素赋值0
                    for n = 1:SIDn(d,a)
                        ARTf(SIDi(n,cd)) = 0;
                    end
                else
                    PDe = 0;  %投影值/投影值差值
                    wt2 = 0;  %权重平方和
                    for n = 1:SIDn(d,a)
                        PDe = PDe+ARTf(SIDi(n,cd))*SIDw(n,cd);  %射线投影值
                        wt2 = wt2+SIDw(n,cd)^2;
                    end
                    PDe = RF*(sino(d,a)-PDe)/wt2;  %射线投影值差值
                    for n = 1:SIDn(d,a)  %更新图像
                        ARTf(SIDi(n,cd)) = ARTf(SIDi(n,cd))+SIDw(n,cd)*PDe;
                    end
                end
            end
        end
    end
    ARTf = min(max(ARTf,0),1);  %约束<0和>1的值
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