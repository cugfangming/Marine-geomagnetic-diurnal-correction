function hy=mt1dtmfh(freq,dz,sig,hf)
%有限元法TM模式大地电磁正演模拟问题
%差分法。
%系统方程采用LU分解法求解。
%注:
%1，在该代码中，假设EXP（-JWT）与时间相关。
%FREQ是真实的频率。
%DZ是一个矢量，包含每层的厚度，排除空气层。
%SIG是一个矢量，包含每层的电导率，不包括空气层。
%检查输入参数的可接受数量
if (nargin < 4)
   error('Not enough input arguments.');
end
if (nargin > 4) 
   error('Too many input arguments.');
end
%   Constants.
mu = 4.0E-7 * pi;
omega = 2.0 * pi * freq;
II = sqrt(-1);
%   扩展模型.
nz = length(sig);
sig(nz+1) = sig(nz);
dz(nz+1) = sqrt(2.0/(sig(nz+1)*omega*mu));
%   计算系数.
% diag elements.
for ki = 1:nz
    diagA(ki) = II * omega * mu * (dz(ki)+dz(ki+1)) ...
                - 2.0/(sig(ki)*dz(ki)) - 2.0/(sig(ki+1)*dz(ki+1));
end
% off-diag elements.
for ki = 1:nz-1
    offdiagA(ki) = 2.0/(sig(ki+1)*dz(ki+1));
end
% 系统矩阵.
mtxA = sparse(diag(diagA) + diag(offdiagA,-1) + diag(offdiagA,1));
% 计算右侧.
% 使用边界条件: hy(1)=hf, hy(nz)=0.0;
rhs = zeros(nz,1);
rhs(1) = -hf*2.0/(sig(1)*dz(1));
%   使用LU分解法求解系统方程： mtxA * hy = rhs 
[L,U,P] = lu(mtxA);
hy = U\(L\P*rhs);
hy = [hf;hy];
% 函数结束.