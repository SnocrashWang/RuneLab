function [omega] = FFT(spd)
%% 带窗傅里叶变换
    len = length(spd);
    w = hamming(len);
    Y = fft((spd-1.1775).*w', len*50, 2);
    P2 = abs(Y/len/50);
    P1 = P2(:,1:len*50/2+1);
    P1(:,2:end-1) = 2*P1(:,2:end-1);
    f = 100*(0:(len*50/2))/len/50 * 2 * pi;

    [~,p] = max(P1(1:300));
    omega = f(p);

%     figure;
%     plot(f(1:300), P1(1:300), "red");

end
