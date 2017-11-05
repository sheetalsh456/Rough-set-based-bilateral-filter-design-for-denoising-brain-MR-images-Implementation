function output = trilateral(noisy_image, N, granule_size)
[row, col] = size(noisy_image);
window_size = 1;
thresh = get_thresholds(noisy_image, N);
%disp(thresh);
class = get_class(noisy_image, row, col, thresh);
rem = get_rem(noisy_image, row, col, granule_size, thresh);
%imshow(rem);
domain_filter = get_domain_filter(window_size);
trilateral_filter = get_trilateral_filter(noisy_image, row, col, window_size,domain_filter, class, rem);
output = trilateral_filter;

function output = get_thresholds(noisy_image, N)
thresh = multithresh(noisy_image,N);
output = thresh;

function output = get_class(noisy_image, row, col, thresh)
class = ones(row, col);
cnt = 1;
for t = thresh
    cnt = cnt + 1;
    for i = 1:row
        for j = 1:col
            if noisy_image(i,j) >= t
                class(i,j) = cnt;
            end
        end
    end
end
output = class;

function output = get_rem(noisy_image, row, col, granule_size, thresh)
rem = zeros(row, col);
for t = thresh
    for i = 1:granule_size:row
        for j = 1:granule_size:col
            all0 = 1;
            all1 = 1;
            for x = 0:granule_size-1
                for y = 0:granule_size-1
                    if i+x <= row && j+y <= col
                        %disp(I(i+x,j+y));
                        if noisy_image(i+x,j+y) <= t
                            all1 = 0;
                        end
                        if noisy_image(i+x,j+y) > t
                            all0 = 0;
                        end
                    end
                end
            end
            if all0 == 0 && all1 == 0
                for x = 0:granule_size-1
                    for y = 0:granule_size-1
                        if i+x <= row && j+y <= col
                           rem(i+x, j+y) = 1;
                        end
                    end
                end 
            else
                for x = 0:granule_size-1
                    for y = 0:granule_size-1
                        if i+x <= row && j+y <= col
                           rem(i+x, j+y) = 0;
                        end
                    end
                end 
            end
        end
    end
end
output = rem;

function output = get_domain_filter(window_size)
sigma_s = 3; 
[x, y]=meshgrid(-window_size:window_size,-window_size:window_size);
domain_filter=exp(-(x.^2+y.^2)/(2*sigma_s^2));  
output = domain_filter;

function output = get_trilateral_filter(noisy_image, row, col, window_size,domain_filter, class, rem)
sigma_r = 0.2;
noisy_image = double(noisy_image)/255.0;
output = zeros(row, col);
rho = zeros(row, col);
rhos = 0.20; 
rhom = 0.55; 
rhol = 0.85;
for i = 1:row
    for j = 1:col
        imin=max(i-window_size,1);
        imax=min(i+window_size,row);
        jmin=max(j-window_size,1);
        jmax=min(j+window_size,col);
        I1=noisy_image(imin:imax,jmin:jmax);
        for x = max(i-window_size,1):min(i+window_size,row)
            for y = max(j-window_size,1):min(j+window_size,col)
                if rem(x, y) == 1
                    if class(x, y) == class(i, j)
                        rho(x, y) = rhom;
                    else
                        rho(x, y) = rhos;
                    end
                else
                    if class(x, y) == class(i, j)
                        if rem(i, j) == 1
                            rho(x, y) = rhom;
                        else
                            rho(x, y) = rhol;
                        end
                    else
                        if rem(i, j) == 1
                            rho(x, y) = rhos;
                        else
                            rho(x, y) = rhos;
                        end
                    end
                end
            end
        end
        rho_filter = rho(imin:imax, jmin:jmax);
        range_filter=exp(-double(I1-noisy_image(i,j)).^2/(2*sigma_r^2)); % Range filter        
        %Taking the product of the range and domain filter.The combination is refered to as Bilater Filter
        bilateral_filter=range_filter.*domain_filter((imin:imax)-i+window_size+1,(jmin:jmax)-j+window_size+1);
        trilateral_filter= bilateral_filter.*rho_filter;       
        Fnorm=sum(trilateral_filter(:));
        output(i,j)=sum(sum(trilateral_filter.*double(I1)))/Fnorm; %normalize the output    
    end
end

